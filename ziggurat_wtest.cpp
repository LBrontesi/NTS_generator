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

double tilted_tempered_stable(XorwowState& s, BitPool& h, double lambda ,double alpha ){
    double gamma = std::pow(lambda,alpha)*alpha*(1-alpha);
    double eps = ((2.0 + sqrt(pi/2.0)) * sqrt(2.0*gamma) + 1.0) / pi;
    double phi = (1/pi) * exp(-(gamma*pi*pi)/8) *(2+sqrt(pi/2))*sqrt(gamma*pi);
    double w1 = eps*sqrt(pi/(2*gamma));
    double w2 = 2*phi*sqrt(pi);
    double w3 = eps*pi;
    double b = (1-alpha)/alpha;
    double U;
    double z;
    double Z;
    double BU;
    for(;;){
        for(;;){
            double V = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;
            double W1 = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;
            if (gamma >=1){
                if (V < w1/(w1+w2)){
                    U = abs(ziggurat(s,h))/sqrt(gamma);
                    if (U >= pi) {
                        continue;
                    }

                } else {
                    U = pi*(1-W1*W1);
                    if (U >= pi) {
                        continue;
                    }

                } 

            } else{

                if (V < w3/(w2+w3)){
                    U = pi*W1;
                    if (U >= pi) {
                        continue;
                    }
               

                } else {
                    U = pi*(1-W1*W1);
                    if (U >= pi) {
                        continue;
                    }

                }

            }
            double W = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;
            BU = std::sin(U) /( std::pow(std::sin(alpha * U), alpha) *std::pow(std::sin((1.0 - alpha) * U), 1.0 - alpha));

            double B0 =
                1.0 /
                (
                    std::pow(alpha, alpha) *
                    std::pow(1.0 - alpha, 1.0 - alpha)
                );

            double c = std::sqrt(BU / B0);
            double vi = pow(sqrt(gamma)+alpha*c,(1/alpha));
            z = vi/(vi-pow(sqrt(gamma),(1/alpha)));
            double gss = 0.0;

            // xi * exp(-gamma * U^2 / 2) * 1[U >= 0, gamma >= 1]
            if (gamma >= 1.0 && U >= 0.0) {
                gss += eps * std::exp(-gamma * U * U / 2.0);
            }

            // psi / sqrt(pi - U) * 1[0 < U < pi]
            if (U > 0.0 && U < pi) {
                gss += phi / std::sqrt(pi - U);
            }

            // xi * 1[0 <= U <= pi, gamma < 1]
            if (gamma < 1.0 && U >= 0.0 && U <= pi) {
                gss += eps;
            }

            double rho =pi *std::exp(
                    -std::pow(lambda, alpha) *
                    (1.0 - 1.0 / (c * c))) *gss/(
                    (1.0 + std::sqrt(pi / 2.0))
                    * std::sqrt(gamma) / c+ z);
            if ((U < pi) && (W*rho<=1)){
                Z = W*rho;
                break;
                
            }
        }
        double a = pow(BU, -1.0/(1.0-alpha));
        double m = pow((b*lambda)/a,alpha);
        double h1 = sqrt((m*alpha)/a);
        double a1 = h1 * sqrt(pi/2);
        double a2 = h1;
        double a3 = z/a;
        double s1 = a1+a2+a3;
        double V1 = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;
        double E1 = 0.0;
        double N1 = 0.0;
        double X;
        if (V1 < a1/s1) {
            N1 = ziggurat(s,h);
            X = m - h1*abs(N1);

        }else if (V1 < (a1+a2)/s1){
            X = ((static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0)*(m+h1-m) + m;

        } else {
            E1 = - log((static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0);
            X = m + h1 + a3 * E1;
        }
        double E = - log(Z);
        if (X <= 0.0) {
            continue;
        }

        double gf = a * (X - m)+ lambda * (std::pow(X, -b) - std::pow(m, -b));
        if (X<m) gf -= (N1*N1)/2;
        if (X>m+h1) gf -= E1;

        if((X>=0) && (gf <= E) ){
            return 1/pow(X,b);
        }
    }
    


}

double nts(XorwowState& state,BitPool& pool,double alpha,double lambda,double beta,double mu,double sigma) {
    
    double T = tilted_tempered_stable(state, pool, lambda, alpha);

    double meanT = alpha * std::pow(lambda, alpha - 1.0);

    double T_scaled = T / meanT;
    double Z = ziggurat(state, pool);

    return mu + beta *( T_scaled-1) + sigma * std::sqrt(T_scaled) * Z;
}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cout << "Usage: ./program <seed>\n";
        return 1;
    }

    const int N = 100000000;

    std::vector<double> numbers(N);

    uint64_t seed = std::stoull(argv[1]);
    BitPool h;
    XorwowState g;
    init_xorwow(g,seed);
    compute_wki();

    auto start3 = std::chrono::high_resolution_clock::now();
    for (int i =0;i <N;i++){
        numbers[i]=nts(g,h,0.5,0.5,0,0.0,1.0);
    }
    auto end3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed3 = end3 - start3;

    std::cout << elapsed3.count() << " ms\n";

    long double sum = 0.0L;

    for (int i = 0; i < N; i++) {
        sum += numbers[i];
    }

    long double mean = sum / N;

    std::cout << "Mean :" << mean << "\n";

    long double sum_sq = 0.0L;


    for (int i = 0; i < N; i++) {
        long double diff = numbers[i]- mean;
        sum_sq += diff*diff;
        
    }

    long double var = sum_sq / N;

    std::cout << "Variance: " << var << "\n";

    long double sum_3 = 0.0L;


    for (int i = 0; i < N; i++) {
        long double diff = numbers[i]- mean;
        sum_3 += diff*diff*diff;
        
    }

    long double skew = (sum_3 / N)/(var*sqrt(var));

    std::cout << "Skewness: " << skew << "\n";

    long double sum_4 = 0.0L;


    for (int i = 0; i < N; i++) {
        long double diff = numbers[i]- mean;
        sum_4 += diff*diff*diff*diff;
        
    }

    long double kur = (sum_4 / N)/(var*var);

    std::cout << "Kurtosis: " << kur << "\n";

    uint32_t count1 = 0;
    uint32_t count2 = 0;
    uint32_t count3 = 0;
    uint32_t count4 = 0;

    for (int i = 0; i < N; i++) {
        count1 += (std::abs(numbers[i]) > 1.0);
        count2 += (std::abs(numbers[i]) > 2.0);
        count3 += (std::abs(numbers[i]) > 3.0);
        count4 += (std::abs(numbers[i]) > 4.0);
    }

    double frequency1 = static_cast<double>(count1) / N;
    double frequency2 = static_cast<double>(count2) / N;
    double frequency3 = static_cast<double>(count3) / N;
    double frequency4 = static_cast<double>(count4) / N;

    std::cout << "Frequency above |1|: " << frequency1 << "\n";
    std::cout << "Frequency above |2|: " << frequency2 << "\n";
    std::cout << "Frequency above |3|: " << frequency3 << "\n";
    std::cout << "Frequency above |4|: " << frequency4 << "\n";

    for (int lag : {1, 2, 5, 10, 50, 100}) {

        long double num = 0.0L;

        for (int i = 0; i < N - lag; i++) {
            num += (numbers[i] - mean)
                * (numbers[i + lag] - mean);
        }

        long double rho = num / sum_sq;

        std::cout << "Lag-" << lag
                << " autocorrelation: "
                << rho << "\n";
    }


    return 0;
}