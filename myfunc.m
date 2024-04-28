function [Err, k, xk_now] = myfunc()
        % Initialize
    mu = .5;
    rho0 = .3;
    Termination = 10e-3;
    Error = Inf;
    k = 1;
    x0 = .3;
    x1 = .7;
    rho = rho0;
    t = 0;
    z = 1;
    xi = 1.7;
    x_star = 1/sqrt(exp(1));
    % This is A Split Line

    xk_now = x1;
    xk_pre = x0;
    while Error >= Termination
        
        alpha = 1/(k+1);
        theta = 1/(2*k+3);
        lambda = 1/(k*sqrt(k));

        % This is A Split Line
        if t ~= z
            x_innertial = inv_exp(xk_now, xk_pre);
            w = my_exp(-alpha, x_innertial, xk_now);

            w_innertial = inv_exp(w, xk_pre);
            z = my_exp(-theta, w_innertial, w);

            t = (z/(exp(xi*rho)))^(1/(1+xi*rho));

        end

        % if t == z
        %     break;
        % end

        trans = my_phi(z) - my_psi(t);
        xk_next = my_exp(rho, trans, t);

        % update rho

        if my_phi(z) - my_phi(t) ~=0
            rho = min(mu*dist(z, t)*dist(t, xk_next)/abs(my_phi(z) - my_phi(t)), rho + lambda);
        else
            rho = rho + lambda;
        end

        Error = dist(xk_next, xk_now);
        Err(k) = Error;
        % Err = Err';
        tt = xk_now;
        xk_now = xk_next;
        xk_pre = tt;
        k = k+1;
  
    end
end