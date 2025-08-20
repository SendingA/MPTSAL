function x0 = solve_triangle(theta1,theta2,phi1,phi2,rTOF)
    theta1 = deg2rad(theta1);
    theta2 = deg2rad(theta2);
    phi1 = deg2rad(phi1);
    phi2 = deg2rad(phi2);
    %假定反射点的坐标是(0,0)
    B = abs(theta1 - theta2); 
    A = abs(phi1 - phi2);
    kb = sin(B)/sin(A); %b = ka
    C = pi - A - B;
    kc = sin(C)/sin(A);
    kd = (kb+kc-1);
    l = 0;r = 1e9;mid = (l+r)/2;
    c = physconst("Lightspeed");
    for i=0:300
        mid = (l+r) /2;
        if(kd * mid < rTOF)
            l = mid;
        else
            r = mid;
        end
    end
    a_len = mid;
    c_len = kc * mid;
    x_sending = a_len * sin(theta1);
    y_sending = a_len * cos(theta1);
    x_reflection = c_len * sin(theta2);
    y_reflection = c_len * cos(theta2);
    x0 = [x_sending,y_sending;x_reflection,y_reflection];
end