function rho = H_P_atomsphere(r_r,r,Sun_rr_Hel)
    
global H_P
    h=r-6.378136300e6;%ae = 6.378136300e6

    if h<110000||h>2000000
        rho = 0;
    else
        count = 1;
        while h>H_P(count,1)
            count = count + 1;
        end
        h_min = H_P(count-1,1);
        h_max = H_P(count,1);

        H_min = (h_min - h_max) / log(H_P(count-1,2)/H_P(count,2));
        H_max = (h_min - h_max) / log(H_P(count,3)/H_P(count,3));

        rho_min = H_P(count,2)*exp((h_min-h)/H_min);
        rho_max = H_P(count,3)*exp((h_min-h)/H_max);

        sindelta=r_r(9)/r;
        cosdelta=sqrt(r_r(7)^2+r_r(8)^2)/r;
        alpha=atan(r_r(8)/r_r(7));
        sindelta_sun=Sun_rr_Hel(3)/norm(Sun_rr_Hel);
        cosdelta_sun=sqrt(Sun_rr_Hel(1)^2+Sun_rr_Hel(2)^2)/norm(Sun_rr_Hel);
        alpha_sun=atan(Sun_rr_Hel(2)/Sun_rr_Hel(1));
        cospsi = sindelta*sindelta_sun + cosdelta*cosdelta_sun*cos(alpha-alpha_sun-deg2rad(30));
        cos_n_psi_2 = ((1+cospsi)/2)^(7/2);

        rho = rho_min + (rho_max - rho_min)*cos_n_psi_2;
    end
    
end

