#include <iostream>
#include <vector>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif


// extern "C"로 선언하여 pRNE_debug를 외부에서 호출할 수 있도록 설정
extern "C" {
    int pRNE_debug(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
    int pRNE_debug_alloc_mem();
    void pRNE_debug_init_mem(int mem);
    void pRNE_debug_free_mem(int mem);

    extern casadi_int pRNE_debug_SZ_ARG;  // 16
    extern casadi_int pRNE_debug_SZ_RES;  // 9
    extern casadi_int pRNE_debug_SZ_IW;   // 0
    extern casadi_int pRNE_debug_SZ_W;    // 2107
}

void pRNE(double* thetalist_in, double* dthetalist_in, double* ddthetalist_in,
          double* thetalist_ref_in, double* dthetalist_ref_in, double* ddthetalist_ref_in,
          double* g_in, double* taulist_out) {

    casadi_int sz_arg = pRNE_debug_SZ_ARG;    
    casadi_int sz_res = pRNE_debug_SZ_RES;
    casadi_int sz_iw = pRNE_debug_SZ_IW;
    casadi_int sz_w = pRNE_debug_SZ_W;

    // 입력 및 출력 배열 설정
    const casadi_real* arg[sz_arg];
    casadi_real* res[sz_res];
    casadi_int iw[sz_iw];
    casadi_real w[sz_w];

    // 입력 배열 초기화
    casadi_real thetalist[7];
    casadi_real dthetalist[7];
    casadi_real ddthetalist[7];
    casadi_real thetalist_ref[7];
    casadi_real dthetalist_ref[7];
    casadi_real ddthetalist_ref[7];
    casadi_real g[3];
    casadi_real Ftip[6] = {0, 0, 0, 0, 0, 0};
    casadi_real Mlist[16 * 8] = { 1.0000, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 1.0000, 0, -0.0002, -0.0431, 0.2097, 1.0000, 0, 0, -1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, 0, 0.0002, -0.1452, 0.3024, 1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, -1.0000, 0, 0, 0, -0.2494, 0.0001, -0.0657, 1.0000, 0, 0, -1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, 0, 0.0002, 0.1139, 0.1561, 1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, -1.0000, 0, 0, 0, -0.1750, -0.0005, 0.0663, 1.0000, 0, 0, -1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, 0, 0.0000, -0.1078, 0.1050, 1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, -1.0000, 0, 0, 0, -0.1008, 0.0003, 0.0047, 1.0000, 1.0000, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 1.0000, 0, -0.0001, 0.0005, 0.0292, 1.0000};;  // 변환 행렬 리스트
    casadi_real Glist[36 * 7] = {
        11.4444, 0, 0, 0, 0, 0, 0, 11.4444, 0, 0, 0, 0, 0, 0, 11.4444, 0, 0, 0, 0, 0, 0, 0.3507, 0.0001, -0.0004, 0, 0, 0, 0.0001, 0.3048, -0.1098, 0, 0, 0, -0.0004, -0.1098, 0.0600, 4.4331, 0, 0, 0, 0, 0, 0, 4.4331, 0, 0, 0, 0, 0, 0, 4.4331, 0, 0, 0, 0, 0, 0, 0.0374, 0.0000, -0.0791, 0, 0, 0, 0.0000, 0.3212, -0.0000, 0, 0, 0, -0.0791, -0.0000, 0.2940, 2.8705, 0, 0, 0, 0, 0, 0, 2.8705, 0, 0, 0, 0, 0, 0, 2.8705, 0, 0, 0, 0, 0, 0, 0.0410, 0.0000, 0.0000, 0, 0, 0, 0.0000, 0.0211, 0.0182, 0, 0, 0, 0.0000, 0.0182, 0.0228, 2.6821, 0, 0, 0, 0, 0, 0, 2.6821, 0, 0, 0, 0, 0, 0, 2.6821, 0, 0, 0, 0, 0, 0, 0.0162, -0.0001, 0.0334, 0, 0, 0, -0.0001, 0.1136, -0.0000, 0, 0, 0, 0.0334, -0.0000, 0.1002, 2.1299, 0, 0, 0, 0, 0, 0, 2.1299, 0, 0, 0, 0, 0, 0, 2.1299, 0, 0, 0, 0, 0, 0, 0.0280, 0.0000, -0.0000, 0, 0, 0, 0.0000, 0.0144, -0.0127, 0, 0, 0, -0.0000, -0.0127, 0.0150, 2.2241, 0, 0, 0, 0, 0, 0, 2.2241, 0, 0, 0, 0, 0, 0, 2.2241, 0, 0, 0, 0, 0, 0, 0.0111, 0.0001, -0.0148, 0, 0, 0, 0.0001, 0.0370, -0.0000, 0, 0, 0, -0.0148, -0.0000, 0.0275, 0.3825, 0, 0, 0, 0, 0, 0, 0.3825, 0, 0, 0, 0, 0, 0, 0.3825, 0, 0, 0, 0, 0, 0, 0.0008, 0, 0.0000, 0, 0, 0, 0, 0.0008, -0.0000, 0, 0, 0, 0.0000, -0.0000, 0.0006};  // 관성 행렬 리스트
    casadi_real Slist[6 * 7] = {
        0, 0, 0, 0, 0, 1.0000, 0.3000, 0, 0, 0, -1.0000, 0, -0.1940, 0, 0, 0, 0, 1.0000, 0.7495, 0, 0, 0.0000, -1.0000, 0, -0.0040, 0, 0, 0, -0.0000, 1.0000, 1.0995, 0.0000, 0, 0.0000, -1.0000, -0.0000, -0.1870, 0, 0, 0, -0.0000, 1.0000};  // 스크류 행렬

    // 입력값 배열 할당
    for (int i = 0; i < 7; ++i) {
        thetalist[i] = thetalist_in[i];
        dthetalist[i] = dthetalist_in[i];
        ddthetalist[i] = ddthetalist_in[i];
        thetalist_ref[i] = thetalist_ref_in[i];
        dthetalist_ref[i] = dthetalist_ref_in[i];
        ddthetalist_ref[i] = ddthetalist_ref_in[i];
    }
    g[0] = g_in[0];
    g[1] = g_in[1];
    g[2] = g_in[2];

    // 포인터 설정
    arg[0] = thetalist;
    arg[1] = dthetalist;
    arg[2] = ddthetalist;
    arg[3] = thetalist_ref;
    arg[4] = dthetalist_ref;
    arg[5] = ddthetalist_ref;
    arg[6] = g;
    arg[7] = Ftip;
    arg[8] = Mlist;
    arg[9] = Glist;
    arg[10] = Slist;

    // 결과값 배열 설정
    casadi_real taulist[7];
    res[0] = taulist;

    // 메모리 할당 및 초기화
    int mem = pRNE_debug_alloc_mem();
    pRNE_debug_init_mem(mem);

    // 함수 호출
    int status = pRNE_debug(arg, res, iw, w, mem);
    if (status != 0) {
        std::cerr << "pRNE_debug 함수 실행에 실패했습니다." << std::endl;
        pRNE_debug_free_mem(mem);
        return;
    }

    // 결과 복사
    for (int i = 0; i < 7; ++i) {
        taulist_out[i] = taulist[i];
    }

    // 메모리 해제
    pRNE_debug_free_mem(mem);
}
