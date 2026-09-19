#include <iostream>
#include <cstdint>
#include <cmath>

/* define costants for the ziggurat algorithm with 256 layers
r is the last right edge, a quantile
v is the area of each rectangle which must be the same*/
#define r 3.6541528853610088
#define v 0.00492867323399

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

/* computes the w, k and f according otthe paper*/
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


double ziggurat(XorwowState& s,BitPool& h){
    for(;;){
        uint32_t j = xorwow(s);
        double sign = random_bit(s,h) ? 1.0 : -1.0;
        uint8_t i = (j&255);
        double x = j*w[i];
        if (j<k[i]){
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

            return sign *(r + xt);
        }
        double u = (static_cast<double>(xorwow(s)) + 0.5) / 4294967296.0;
        if((u*(f[i-1]-f[i]))<(exp(-0.5*x*x)-f[i])){
            return sign*x;
        } 
        

    }
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
    auto start = std::chrono::high_resolution_clock::now();
    for (int i =0;i<N;i++){
        numbers[i] = ziggurat(g,h);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << elapsed.count() << " ms\n";

    long double sum = 0.0L;

    for (int i = 0; i < N; i++) {
        sum += numbers[i];
    }

    long double mean = sum / N;

    std::cout << "Mean: " << mean << "\n";

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