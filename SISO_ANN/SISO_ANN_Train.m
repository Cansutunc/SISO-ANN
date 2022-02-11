function [X,y_m_value]=SISO_ANN_Train(s,t,y_real)
    X=-2+4*rand(3*s+1,1); %Weight matrix with random values is created
    syms x [3*s+1 1] %Symbolic expressions created
    h=(exp(x1)-exp(-x1))/(exp(x1)+exp(-x1)); %Activation Function
    y_model=vpa(zeros(size(y_real,1),1)); %Empty Matrix created for Output of Model
    J=vpa(zeros(size(y_real,1),3*s+1)); %Empty Matrix created for Jacobian Matrix
    o_max=30; %Maximum Iteration Number for ANN
    o=1; %Iteration Counter
    
    %Calculation of Model Outputs and Jacobian
    for i=1:size(t,1)
        %Calculation of Model Outputs as a Symbolic Function
        sum=0;
        for j=(2*s+1):3*s
            mult=x(j)*subs(h,x1,x(j-2*s)*t(i)+x(j-s));
            sum=mult+sum;
        end
        y_model(i,1)=sum+x(3*s+1);
       
        %Calculation of Jacobian Matrix
        for j=1:s
            h_v=subs(h,x1,x(j)*t(i)+x(s+j));
            J(i,j)=-(x(j+2*s)*t(i)*(1-h_v^2));
        end
        for j=s+1:2*s
            h_v=subs(h,x1,x(j-s)*t(i)+x(j));
            J(i,j)=-(x(j+s)*(1-h_v^2));
        end
        for j=2*s+1:3*s
            J(i,j)=-(subs(h,x1,x(j-2*s)*t(i)+x(j-s)));
        end
        for j=3*s+1
            J(i,j)=-1;
        end
     
    end
    E=vpa(y_real-y_model); %Error Function
    F=vpa(E.'*E); %Function to Optimize to find best Weight Values
    while 1
        n_max=100; %Number of Maximum Iteration for LM.
        mu=1; %Mu (Between mu_min and mu_max) 
        mu_scale=10; %Mu scale
        mu_min=0.00001; %Mu minimum
        mu_max=10000; %Mu maximum
        E1=10^-4; %Terminate Condition 1 for LM
        E2=10^-4; %Terminate Condition 2 for LM
        E3=10^-4; %Terminate Condition 3 for LM
        X_old=X; %Old Weight Matrix hold in this variable to see progress of the modelling
        [X]=Levenberg_M(F,X,E,J,mu,mu_scale,mu_min,mu_max,n_max,E1,E2,E3); %LM optimization called
        fprintf('Number of Iteration of ANN is %d  right now\n',o); %Number of Iteration of ANN Displayed to User
        norm_E=double(norm(subs(E,x,X))); %Norm of Error Matrix Calculated for Ending Condition
        delta_X=norm(X-X_old); %Differences Between the Old and New Weight Matrixes Calculated for Ending Condition
        f_MSE=double(subs(F,x,X)); %Value of f Function Calculated for Ending Condition
        training_error(o)=log10(norm_E); %Training Error Calculated for Ending Condition
        if norm_E<0.035 && delta_X<5*10^-4 || f_MSE<10^-4|| training_error(o)<-3 || o>o_max
            fprintf('\n<SISO_ANN_Train> Number of Iterations of SISO ANN Training is %d\n',o);
            fprintf('<SISO_ANN_Train> Training Error is %.4f\n',norm_E);
            fprintf('<SISO_ANN_Train> Change in ANN Weights is %.4f\n',delta_X);
            y_m_value=double(subs(y_model,x,X));
            figure(1);
            subplot(2,1,1);
            plot(t,y_real,'-*');
            hold on
            plot(t,y_m_value,'-o');
            grid on
            title('Train Data Graph');
            xlabel('Input');
            ylabel('Output');
            legend('Real Output of Training Data','Model Output of Train Data')
            
            subplot(2,1,2);
            plot(1:1:o,training_error,'b');
            grid on
            title('Train Error');
            xlabel('Iteration');
            ylabel('log(MSE)');
            break;
        end
        o=o+1;
    end
end