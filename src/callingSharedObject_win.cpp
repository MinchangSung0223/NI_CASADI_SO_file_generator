#include <iostream>
#include <Eigen/Dense>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// CasADi에서 생성된 함수 포인터 정의
typedef double casadi_real;
typedef int casadi_int;

int main()
{
#ifdef _WIN32
    // Windows 환경에서의 라이브러리 로드
    HMODULE handle = LoadLibrary(TEXT("pRNE.dll")); // 공유 라이브러리 경로 설정
    if (!handle)
    {
        std::cerr << "Failed to load library." << std::endl;
        return EXIT_FAILURE;
    }

    // GravityForces 함수 포인터 가져오기
    typedef int (*GravityForces_t)(const casadi_real **arg, casadi_real **res, casadi_int *iw, casadi_real *w, int mem);
    GravityForces_t GravityForces = (GravityForces_t)GetProcAddress(handle, "GravityForces");
    if (!GravityForces)
    {
        std::cerr << "Failed to load GravityForces function." << std::endl;
        FreeLibrary(handle);
        return EXIT_FAILURE;
    }

    // GravityForces_work 함수 포인터 가져오기
    typedef int (*GravityForces_work_t)(casadi_int *sz_arg, casadi_int *sz_res, casadi_int *sz_iw, casadi_int *sz_w);
    GravityForces_work_t GravityForces_work = (GravityForces_work_t)GetProcAddress(handle, "GravityForces_work");

#else
    // Linux 환경에서의 라이브러리 로드
    void *handle = dlopen("./libmylib_ubuntu.so", RTLD_LAZY); // 공유 라이브러리 경로 설정
    if (!handle)
    {
        std::cerr << "Failed to load library: " << dlerror() << std::endl;
        return EXIT_FAILURE;
    }

    // GravityForces 함수 포인터 가져오기
    typedef int (*GravityForces_t)(const casadi_real **arg, casadi_real **res, casadi_int *iw, casadi_real *w, int mem);
    GravityForces_t GravityForces = (GravityForces_t)dlsym(handle, "GravityForces");
    if (!GravityForces)
    {
        std::cerr << "Failed to load GravityForces function: " << dlerror() << std::endl;
        dlclose(handle);
        return EXIT_FAILURE;
    }

    // GravityForces_work 함수 포인터 가져오기
    typedef int (*GravityForces_work_t)(casadi_int *sz_arg, casadi_int *sz_res, casadi_int *sz_iw, casadi_int *sz_w);
    GravityForces_work_t GravityForces_work = (GravityForces_work_t)dlsym(handle, "GravityForces_work");
#endif

    // Work size 확인 및 메모리 할당
    casadi_int sz_arg, sz_res, sz_iw, sz_w;
    GravityForces_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

    casadi_real **arg = (casadi_real **)malloc(sz_arg * sizeof(casadi_real *));
    casadi_real **res = (casadi_real **)malloc(sz_res * sizeof(casadi_real *));
    casadi_real *w = (casadi_real *)malloc(sz_w * sizeof(casadi_real));

    // thetalist 값 설정 (예: 7개 값)
    double thetalist[7] = {0.5, 0.7, 0.5, 0.2, 0.1, 0.3, 0.6};
    double g[3] = {0.0, 0.0, -9.8}; // 중력 벡터

    // 주어진 Mlist_array 값
     double Mlist_array[128] = {
        1.0000, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 1.0000, 0, -0.0002, -0.0431, 0.2097, 1.0000, 0, 0, -1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, 0, 0.0002, -0.1452, 0.3024, 1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, -1.0000, 0, 0, 0, -0.2494, 0.0001, -0.0657, 1.0000, 0, 0, -1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, 0, 0.0002, 0.1139, 0.1561, 1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, -1.0000, 0, 0, 0, -0.1750, -0.0005, 0.0663, 1.0000, 0, 0, -1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, 0, 0.0000, -0.1078, 0.1050, 1.0000, 0, 1.0000, 0, 0, 0, 0, -1.0000, 0, -1.0000, 0, 0, 0, -0.1008, 0.0003, 0.0047, 1.0000, 1.0000, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 1.0000, 0, -0.0001, 0.0005, 0.0292, 1.0000};

    // 주어진 Glist_array 값
    double Glist_array[252] = {
        11.4444, 0, 0, 0, 0, 0, 0, 11.4444, 0, 0, 0, 0, 0, 0, 11.4444, 0, 0, 0, 0, 0, 0, 0.3507, 0.0001, -0.0004, 0, 0, 0, 0.0001, 0.3048, -0.1098, 0, 0, 0, -0.0004, -0.1098, 0.0600, 4.4331, 0, 0, 0, 0, 0, 0, 4.4331, 0, 0, 0, 0, 0, 0, 4.4331, 0, 0, 0, 0, 0, 0, 0.0374, 0.0000, -0.0791, 0, 0, 0, 0.0000, 0.3212, -0.0000, 0, 0, 0, -0.0791, -0.0000, 0.2940, 2.8705, 0, 0, 0, 0, 0, 0, 2.8705, 0, 0, 0, 0, 0, 0, 2.8705, 0, 0, 0, 0, 0, 0, 0.0410, 0.0000, 0.0000, 0, 0, 0, 0.0000, 0.0211, 0.0182, 0, 0, 0, 0.0000, 0.0182, 0.0228, 2.6821, 0, 0, 0, 0, 0, 0, 2.6821, 0, 0, 0, 0, 0, 0, 2.6821, 0, 0, 0, 0, 0, 0, 0.0162, -0.0001, 0.0334, 0, 0, 0, -0.0001, 0.1136, -0.0000, 0, 0, 0, 0.0334, -0.0000, 0.1002, 2.1299, 0, 0, 0, 0, 0, 0, 2.1299, 0, 0, 0, 0, 0, 0, 2.1299, 0, 0, 0, 0, 0, 0, 0.0280, 0.0000, -0.0000, 0, 0, 0, 0.0000, 0.0144, -0.0127, 0, 0, 0, -0.0000, -0.0127, 0.0150, 2.2241, 0, 0, 0, 0, 0, 0, 2.2241, 0, 0, 0, 0, 0, 0, 2.2241, 0, 0, 0, 0, 0, 0, 0.0111, 0.0001, -0.0148, 0, 0, 0, 0.0001, 0.0370, -0.0000, 0, 0, 0, -0.0148, -0.0000, 0.0275, 0.3825, 0, 0, 0, 0, 0, 0, 0.3825, 0, 0, 0, 0, 0, 0, 0.3825, 0, 0, 0, 0, 0, 0, 0.0008, 0, 0.0000, 0, 0, 0, 0, 0.0008, -0.0000, 0, 0, 0, 0.0000, -0.0000, 0.0006};

    // 주어진 Slist_array 값
    double Slist_array[42] = {
        0, 0, 0, 0, 0, 1.0000, 0.3000, 0, 0, 0, -1.0000, 0, -0.1940, 0, 0, 0, 0, 1.0000, 0.7495, 0, 0, 0.0000, -1.0000, 0, -0.0040, 0, 0, 0, -0.0000, 1.0000, 1.0995, 0.0000, 0, 0.0000, -1.0000, -0.0000, -0.1870, 0, 0, 0, -0.0000, 1.0000};
    // CasADi에 전달할 입력 배열 설정
    arg[0] = thetalist;   // thetalist (7)
    arg[1] = g;           // 중력 벡터 (3)
    arg[2] = Mlist_array; // Mlist (16x8 -> 128)
    arg[3] = Glist_array; // Glist (36x7 -> 252)
    arg[4] = Slist_array; // Slist (6x7 -> 42)

    // 출력 배열 설정 (7개 결과 값)
    res[0] = (casadi_real *)malloc(7 * sizeof(casadi_real));

    // GravityForces 함수 호출
    if (GravityForces((const casadi_real **)arg, res, nullptr, w, 0) != 0)
    {
        std::cerr << "GravityForces function failed." << std::endl;
#ifdef _WIN32
        FreeLibrary(handle);
#else
        dlclose(handle);
#endif
        return EXIT_FAILURE;
    }

    // 결과 출력
    std::cout << "Gravity Forces Result: " << std::endl;
    for (int i = 0; i < 7; ++i)
    {
        std::cout << "res[0][" << i << "] = " << res[0][i] << std::endl;
    }

    // 메모리 해제
    free(res[0]);
    free(arg);
    free(res);
    free(w);

    // 동적 라이브러리 닫기
#ifdef _WIN32
    FreeLibrary(handle);
#else
    dlclose(handle);
#endif

    return 0;
}
