function [X]=LM(f,X,E,J,mu,mu_scale,mu_min,mu_max,n_max,E1,E2,E3)
    syms x [1 length(X)];
    grad_F=2.*J.'*E;
    X_old=0;
    c=size(J,2);
    i=1;
    while 1
        f_value=double(subs(f,x.',X)); 
        J_value=double(subs(J,x.',X)); 
        E_value=double(subs(E,x.',X)); 
        while 1
            zk=-(inv(J_value.'*J_value+mu.*eye(c)))*J_value.'*E_value;
            f_xkzk=double(subs(f,x.',X+zk));
            if f_xkzk<f_value
                pk=zk;
                f_sk=subs(f,x.',X+x1.*pk);
                [sk]=golden_section(f_sk,50,-50,10^-4);
                X_old=X;
                X=double(X+sk.*pk);
                mu=mu/mu_scale;
                break
            else
                mu=mu*mu_scale;
                if mu<mu_max && mu>mu_min
                    continue; 
                else
                    break; 
                end
            end
        end
        f1=subs(f,x.',X_old);
        f2=subs(f,x.',X);
        normgrad=norm(subs(grad_F,x.',X));
        error=abs(f2-f1);
        delta_x=norm(X-X_old);
        if i>n_max
            fprintf('<Levenberg_Marquardt> The number of iterations exceeded the maximum value.\n');
            break
        elseif error<E1
            fprintf('<Levenberg_Marquardt> Condition 2 is established in %d iterations.\n',i);
            break
        elseif delta_x<E2
            fprintf('<Levenberg_Marquardt> Condition 3 is established in %d iterations.\n',i);
            break
        elseif normgrad<E3
            fprintf('<Levenberg_Marquardt> Condition 4 is established in %d iterations.\n',i);
            break
        elseif mu<mu_min || mu>mu_max
            fprintf('<Levenberg_Marquardt> Condition 5 is established in %d iterations.\n',i);
            break;
        end
        i=i+1;
    end
end