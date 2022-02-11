function[xk]=golden_section(f,xh,xl,xend)
    syms x1;
    T=0.38197;
    tol=xend/(xh-xl);
    ni=-2.078*log(tol);
      
    a1=xl+T*(xh-xl);
    f1=subs(f,x1,a1);
    a2=xh-T*(xh-xl);
    f2=subs(f,x1,a2);

    i=1;
    while 1
        if i<ni
            if f1>f2
                xl=a1;
                a1=a2;
                f1=f2;
                a2=xh-T*(xh-xl);
                f2=subs(f,x1,a2); 
            elseif f2>f1
                xh=a2;
                a2=a1;
                f2=f1;
                a1=xl+T*(xh-xl);
                f1=subs(f,x1,a1);
            end
        else
            xk=double((a1+a2)/2);
            fprintf('<GOLDEN SECTÝON>The Root is=%.4f\n',xk);
            break;
        end
        i=i+1;
    end
end