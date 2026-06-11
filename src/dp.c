/* 1+1D directed site percolation with two-valued quenched columnar disorder.
 * Original implementation. Site i active at t+1 iff (i or i+1 active at t) && U < p_i.
 * RNG: xoshiro256++ (public-domain algorithm, implemented here from the paper spec).
 * Bernoulli via 32-bit fixed-point threshold (p resolution ~2.3e-10).
 * Usage: dp <mode 0=clean,1=squarefree,2=shuffled> <p> <c> <L> <T> <offset> <seed> <nsamp_out>
 * Output: t rho pairs at log-spaced times.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static uint64_t s[4];
static inline uint64_t rotl(uint64_t x, int k){ return (x<<k)|(x>>(64-k)); }
static inline uint64_t next(void){
    uint64_t r = rotl(s[0]+s[3],23)+s[0];
    uint64_t t = s[1]<<17;
    s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=t; s[3]=rotl(s[3],45);
    return r;
}
static void seed_rng(uint64_t seed){
    /* splitmix64 init, from the reference description */
    uint64_t z = seed;
    for(int i=0;i<4;i++){
        z += 0x9e3779b97f4a7c15ULL;
        uint64_t x = z;
        x = (x ^ (x>>30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x>>27)) * 0x94d049bb133111ebULL;
        s[i] = x ^ (x>>31);
    }
}

int main(int argc, char **argv){
    if(argc < 9){ fprintf(stderr,"args\n"); return 1; }
    int mode = atoi(argv[1]);
    double p = atof(argv[2]);
    double c = atof(argv[3]);
    long L = atol(argv[4]);
    long T = atol(argv[5]);
    uint64_t offset = strtoull(argv[6],0,10);
    uint64_t seed = strtoull(argv[7],0,10);
    int nout = atoi(argv[8]);
    seed_rng(seed);

    /* squarefree mask of [offset, offset+L) */
    unsigned char *good = malloc(L);
    memset(good,1,L);
    if(mode != 0){
        for(uint64_t d=2; d*d < offset+(uint64_t)L; d++){
            uint64_t q = d*d;
            uint64_t start = ((offset + q - 1)/q)*q;
            for(uint64_t m = start; m < offset+(uint64_t)L; m += q) good[m-offset]=0;
        }
        if(mode == 2){ /* Fisher-Yates shuffle of the mask */
            for(long i=L-1;i>0;i--){
                long j = (long)(next() % (uint64_t)(i+1));
                unsigned char tmp=good[i]; good[i]=good[j]; good[j]=tmp;
            }
        }
    }
    uint32_t thr_good = (uint32_t)(p * 4294967296.0);
    uint32_t thr_bad  = (uint32_t)(c * p * 4294967296.0);
    uint32_t *thr = malloc(L*sizeof(uint32_t));
    for(long i=0;i<L;i++) thr[i] = good[i]? thr_good : thr_bad;

    unsigned char *st = malloc(L), *nx = malloc(L);
    memset(st,1,L);
    /* log-spaced sample times */
    long *ts = malloc(nout*sizeof(long)); int nt=0;
    for(int k=0;k<nout;k++){
        long t = (long)llround(pow(10.0, log10((double)T)*k/(nout-1)));
        if(t<1) t=1;
        if(nt==0 || t>ts[nt-1]) ts[nt++]=t;
    }
    int j=0;
    for(long t=1;t<=T;t++){
        long i=0;
        for(;i+2<=L;i+=2){
            uint64_t r = next();
            unsigned char a0 = (st[i]|st[i+1]) & ((uint32_t)r < thr[i]);
            unsigned char a1 = (st[i+1]|st[(i+2)%L]) & ((uint32_t)(r>>32) < thr[i+1]);
            nx[i]=a0; nx[i+1]=a1;
        }
        for(;i<L;i++){
            uint64_t r = next();
            nx[i] = (st[i]|st[(i+1)%L]) & ((uint32_t)r < thr[i]);
        }
        /* fix periodic wrap for last pair element referencing (i+2)%L==0 done above; redo last site cleanly */
        nx[L-1] = (st[L-1]|st[0]) & ((uint32_t)next() < thr[L-1]);
        unsigned char *tmp=st; st=nx; nx=tmp;
        if(j<nt && t==ts[j]){
            long sum=0; for(long k2=0;k2<L;k2++) sum+=st[k2];
            printf("%ld %.10e\n", t, (double)sum/L);
            j++;
            if(sum==0) break;
        }
    }
    return 0;
}
