function [y_model]=SISO_ANN_Test(X,s,test_data,y_real_test)
    syms x1
    h=(exp(x1)-exp(-x1))/(exp(x1)+exp(-x1));
    y_model=zeros(size(test_data,1),1);
    for i=1:size(test_data,1)
        %Calculation of Model Outputs
        sum=0;
        for j=(2*s+1):3*s
            mult=X(j)*subs(h,x1,X(j-2*s)*test_data(i)+X(j-s));
            sum=mult+sum;
        end
        y_model(i,1)=double(sum+X(3*s+1));
    end
    test_Error=norm(y_model-y_real_test);
    fprintf('\n<SISO_ANN_Test> Test Data Eror is %f\n',test_Error);
    figure(2);
    plot(test_data,y_real_test,'-*');
    hold on
    plot(test_data,y_model,'-o');
    grid on
    title('Test Data Graph');
    xlabel('Input');
    ylabel('Output');
    legend('Real Output of Test Data','Model Output of Test Data')
end