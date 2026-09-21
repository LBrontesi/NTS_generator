#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>

/* define costants for the ziggurat algorithm with 256 layers
r is the last right edge, a quantile
v is the area of each rectangle which must be the same*/
#define r 3.6541528853610088
#define v 0.00492867323399
#define pi 3.1415926535897932

double x[256];
double w[256];
uint32_t k[256];
double f[256];

struct XorwowState {
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint32_t a;
    uint32_t b;
    uint32_t d;   
};

struct BitPool {
    uint32_t bits;
    int remaining = 0;
};

struct NTSparams {
    double alpha;
    double one_minus_alpha;
    double lambda;

    double inv_alpha;
    double inv_one_minus_alpha;
    double neg_inv_one_minus_alpha;

    double lambda_alpha;

    double gamma;
    double gamma_sqrt;
    double gamma_power;

    double eps;
    double phi;

    double w1;
    double w2;
    double w3;

    double prob_w1_w2;
    double prob_w2_w3;

    double b;
    double b_lambda;

    double inv_B0;

    double meanT;
    double inv_meanT;

    double sqrt_pi_over_2;
    double rho_gamma_coeff;
    double mix12_coeff;

    bool gamma_ge_one;
};


/* Function needed to generate a 64 bits unsigned integer used to initialize Xorwow  struct 
given xorwow needs 32 bits for a state, then one generation from splitmix64 can be used to initialize 2 states in xorwow
*/
uint64_t splitmix64(uint64_t& state) {
    uint64_t z = (state += 0x9E3779B97F4A7C15ULL);

    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;

    return z ^ (z >> 31);
}

/* function used to initialize the xorwow struct using the splitmix64
the original seed is not modified because seed in init xorwowis passed by value and not reference
so a copy is created and that copy is modified internally by splitmix64 otherwise it would generate alwyas thensame value
*/

void init_xorwow(XorwowState& g, uint64_t seed) {
    uint64_t sm_state = seed;

    uint64_t r1 = splitmix64(sm_state);
    g.x = static_cast<uint32_t>(r1);
    g.y = static_cast<uint32_t>(r1 >> 32);

    uint64_t r2 = splitmix64(sm_state);
    g.z = static_cast<uint32_t>(r2);
    g.a = static_cast<uint32_t>(r2 >> 32);

    uint64_t r3 = splitmix64(sm_state);
    g.b = static_cast<uint32_t>(r3);
    g.d = static_cast<uint32_t>(r3 >> 32);
}

/*logic of xorwow algorithm, right shift and XOR operation to create a temporary value, this value is mixed with 
actual state to create a new state via left shifts and XOR operations, then Weyl costant is updated
and states are moved left, the value returned is the new state plus Weyl costant*/
uint32_t xorwow(XorwowState& s){
    uint32_t t = s.x^(s.x>>2);
    uint32_t new_state = s.b ^(s.b<<4)^t^(t<<1);
    s.d += 362437;
    s.x = s.y;
    s.y = s.z;
    s.z = s.a;
    s.a = s.b;
    s.b = new_state;
    return new_state + s.d;
}

/*given i need a random sign, this function takes as input a XORWOW generation and keeps 1 bit for the sign
the other 31 bits are kept to generate signs later*/
uint32_t random_bit(XorwowState& s, BitPool& pool) {
    if (pool.remaining == 0) {
        pool.bits = xorwow(s);
        pool.remaining = 32;
    }

    uint32_t bit = pool.bits & 1u;
    pool.bits >>= 1;
    pool.remaining--;

    return bit;
}


/*computes the xi recurseviley based on xi = f^-1(v/xi+1 + f(x+1))
*/
void compute_xi(){
    x[255] = r;
    x[0] = 0;
    for(int i = 254;i>0;i--){
        x[i] = sqrt(-2*log(v/x[i+1]+exp(-(x[i+1]*x[i+1])/2)));
    }  
}

/* computes the w, k and f according to the paper baesed on the xi*/
void compute_wki(){
    compute_xi();
    for(int i = 0;i<256;i++){        
        if(i==0){
            k[i] = static_cast<uint32_t>((4294967296.0 * r*exp(-r*r/2))/v);
            w[i] = ((1/4294967296.0)* (v/exp(-r*r/2)));
            f[i] = 1.0;
        } else{
            k[i] = static_cast<uint32_t>(4294967296.0 * (x[i - 1] / x[i]));
            w[i] = ((1/4294967296.0)* x[i]);
            f[i] = exp(-(x[i]*x[i])/2);
        }
    }  
}

/* implement ofthe alogrithm using a initiliazed XORWOW struct to generate a 32 bits unsigned number, takes the first 8
bits for the index of the layer given they are 256 and follow the paper using first fast acceptance, if rejected
goes to the 
*/
double ziggurat(XorwowState& s,BitPool& h){
    for(;;){
        uint32_t j = xorwow(s);
        uint8_t i = (j&255);
        double x = j*w[i];
        if (j<k[i]){
            double sign = random_bit(s,h) ? 1.0 : -1.0;
            return sign*x;
        }
        if (i == 0) {
            double xt, y;
            do {
                double u1 = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;
                double u2 = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;

                xt = -std::log(u1) / r;
                y  = -std::log(u2);

            } while (2.0 * y < xt * xt);
            double sign = random_bit(s,h) ? 1.0 : -1.0;
            return sign *(r + xt);
        }
        double u = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;
        if((u*(f[i-1]-f[i]))<(exp(-0.5*x*x)-f[i])){
            double sign = random_bit(s,h) ? 1.0 : -1.0;
            return sign*x;
        } 
        

    }
}

void initNTSparams(
    NTSparams& p,
    double lambda,
    double alpha
) {
    p.alpha = alpha;
    p.one_minus_alpha = 1.0 - alpha;
    p.lambda = lambda;

    p.inv_alpha = 1.0 / alpha;

    p.inv_one_minus_alpha =
        1.0 / p.one_minus_alpha;

    p.neg_inv_one_minus_alpha =
        -p.inv_one_minus_alpha;

    // lambda^alpha
    p.lambda_alpha =
        std::pow(lambda, alpha);

    p.gamma =
        p.lambda_alpha
        * alpha
        * p.one_minus_alpha;

    p.gamma_sqrt =
        std::sqrt(p.gamma);

    // (sqrt(gamma))^(1/alpha)
    p.gamma_power =
        std::pow(
            p.gamma_sqrt,
            p.inv_alpha
        );

    p.sqrt_pi_over_2 =
        std::sqrt(pi / 2.0);

    p.eps =
        (
            (2.0 + p.sqrt_pi_over_2)
            * std::sqrt(2.0 * p.gamma)
            + 1.0
        )
        / pi;

    p.phi =
        (1.0 / pi)
        * std::exp(
            -(p.gamma * pi * pi) / 8.0
        )
        * (2.0 + p.sqrt_pi_over_2)
        * std::sqrt(p.gamma * pi);

    p.w1 =
        p.eps
        * std::sqrt(
            pi / (2.0 * p.gamma)
        );

    p.w2 =
        2.0
        * p.phi
        * std::sqrt(pi);

    p.w3 =
        p.eps * pi;

    p.prob_w1_w2 =
        p.w1 / (p.w1 + p.w2);

    p.prob_w2_w3 =
        p.w3 / (p.w2 + p.w3);

    p.b =p.one_minus_alpha / alpha;
    p.b_lambda =p.b * lambda;
    p.inv_B0 = std::pow(alpha, alpha) *std::pow(p.one_minus_alpha,p.one_minus_alpha);
    p.meanT =alpha * std::pow(lambda,alpha - 1.0);
    p.inv_meanT =1.0 / p.meanT;
    p.rho_gamma_coeff =(1.0 + p.sqrt_pi_over_2)* p.gamma_sqrt;
    p.mix12_coeff = p.sqrt_pi_over_2 + 1.0;
    p.gamma_ge_one = p.gamma >= 1.0;
}


inline double uniform01(XorwowState& s) {
    constexpr double INV_2_32 =1.0 / 4294967296.0;
    return(static_cast<double>(xorwow(s)) + 0.5)* INV_2_32;
}

double tilted_tempered_stable( XorwowState& s,BitPool& h,const NTSparams& p) {
    double U;
    double z;
    double Z;
    double BU;
    for (;;) {

        // --------------------------
        // First rejection stage
        // --------------------------

        for (;;) {

            double V  = uniform01(s);
            double W1 = uniform01(s);

            if (p.gamma_ge_one) {

                if (V < p.prob_w1_w2) {
                    U =
                        std::abs(ziggurat(s, h))
                        / p.gamma_sqrt;
                }
                else {
                    U =
                        pi * (1.0 - W1 * W1);
                }

            }
            else {

                if (V < p.prob_w2_w3) {
                    U = pi * W1;
                }
                else {
                    U =
                        pi * (1.0 - W1 * W1);
                }
            }

            if (U >= pi)
                continue;

            double W = uniform01(s);

            BU =
                std::sin(U)
                /
                (
                    std::pow(
                        std::sin(p.alpha * U),
                        p.alpha
                    )
                    *
                    std::pow(
                        std::sin(
                            p.one_minus_alpha * U
                        ),
                        p.one_minus_alpha
                    )
                );

            // division removed
            double c =
                std::sqrt(
                    BU * p.inv_B0
                );

            double vi =
                std::pow(
                    p.gamma_sqrt
                    + p.alpha * c,
                    p.inv_alpha
                );

            // pow(gamma_sqrt,1/alpha)
            // no longer calculated here
            z =
                vi /
                (vi - p.gamma_power);

            double gss;

            if (p.gamma_ge_one) {

                gss =
                    p.eps
                    * std::exp(
                        -p.gamma
                        * U * U
                        / 2.0
                    );

            }
            else {

                gss = p.eps;
            }

            if (U > 0.0) {
                gss +=
                    p.phi
                    / std::sqrt(pi - U);
            }

            double rho =pi*std::exp(-p.lambda_alpha*(1.0- 1.0 / (c * c)) )*gss/(p.rho_gamma_coeff / c + z);

            Z = W * rho;

            if (Z <= 1.0)
                break;
        }

        double a =
            std::pow(
                BU,
                p.neg_inv_one_minus_alpha
            );

        double m =
            std::pow(
                p.b_lambda / a,
                p.alpha
            );

        double h1 =
            std::sqrt(
                (m * p.alpha) / a
            );

        double a1 =
            h1 * p.sqrt_pi_over_2;

        double a3 =
            z / a;

        // a1 + a2 + a3
        // a2 = h1
        double s1 =h1 * p.mix12_coeff+ a3;
        double V1 = uniform01(s);
        double E1 = 0.0;
        double N1 = 0.0;
        double X;

        if (V1 < a1 / s1) {
            N1 = ziggurat(s, h);
            X = m- h1 * std::abs(N1);
        }
        else if (V1 <(h1 * p.mix12_coeff) / s1) {
            X = m + h1 * uniform01(s);
        }
        else {
            E1 =-std::log(uniform01(s));
            X = m+ h1+ a3 * E1;
        }
        if (X <= 0.0)
            continue;

        double E = -std::log(Z);
        double X_neg_b =std::pow(X, -p.b);
        double m_neg_b =std::pow(m, -p.b);
        double gf = a * (X - m)+p.lambda    * (X_neg_b - m_neg_b);

        if (X < m)
            gf -= (N1 * N1) / 2.0;

        if (X > m + h1)
            gf -= E1;

        if (gf <= E) {
            return X_neg_b;
        }
    }
}

double nts( XorwowState& state,BitPool& pool,const NTSparams& p,double beta,double mu,double sigma) {
    double T = tilted_tempered_stable(state,pool,p);

    double T_scaled = T * p.inv_meanT;

    double Z = ziggurat(state,pool);

    return mu+ beta * (T_scaled - 1.0)+ sigma* std::sqrt(T_scaled)* Z;
}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cout << "Usage: ./program <seed>\n";
        return 1;
    }

    const int N = 10000000;

    std::vector<double> numbers(N);

    uint64_t seed = std::stoull(argv[1]);
    BitPool h;
    XorwowState g;
    NTSparams p;
    initNTSparams(p,0.5,0.5);
    init_xorwow(g,seed);
    compute_wki();

    auto start3 = std::chrono::high_resolution_clock::now();
    for (int i =0;i <N;i++){
        numbers[i]=nts(g,h,p,0.5,0.0,1.0);
    }
    auto end3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed3 = end3 - start3;

    std::cout << "Number of RVs generated " << N << " in " << elapsed3.count() << " ms" << "\n";
}