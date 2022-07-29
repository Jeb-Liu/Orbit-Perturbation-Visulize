function sigma_moon_start=moonJ2000(JD_T)
    a = 384747981;
    e = 0.054879905;
    ii = deg2rad( 5.12983501671 );
    W = deg2rad( 125.044556 - 1934.136185*JD_T + 0.0020767*JD_T^2 );
    G = deg2rad( 83.3532417 + 4069.013711*JD_T - 0.0103236*JD_T^2 );
    L = deg2rad( 218.3166556 + 481267.8813425*JD_T - 0.00132972*JD_T^2 );
    
    w = G - W;
    l = L - G;

    sigma_moon_start = [a e ii W w l];
end
