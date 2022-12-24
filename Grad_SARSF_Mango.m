function x=Grad_SARSF_Mango(fun,gfun, x0, W, X, Y, epsilon)
%����: �ݶȷ������Լ���Ż�����:  min f(x)
%����: fun, gfun�ֱ���Ŀ�꺯�������ݶ�, x0�ǳ�ʼ��,
%      epsilonΪ�������
%���: k�ǵ�������, x, val�ֱ��ǽ������ŵ������ֵ
maxk=1000;  %����������
beta=0.1;  sigma=0.4;
k=0;
while(k<maxk)
    gk=feval(gfun,W, X, Y, x0);  %�����ݶ�
    old_ln=feval(fun, W, X, Y, x0);
    dk=-gk;  %������������
    old_norm=norm(gk);
    if(old_norm<epsilon), break; end  %������ֹ׼��
    m=0;  mk=0;
    while(m<20)  %��Armijo�����󲽳�
        temp = x0+beta^m*dk;
        if temp(length(x0)-1)>1 || temp(length(x0)-1)<0 || temp(length(x0))<0
            if m==0
                x = x0;
                return;
            else
                mk=m-1; break;
            end
        end
        
        if(feval(fun,W, X, Y, temp)<=feval(fun,W, X, Y, x0)+sigma*beta^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    x0=x0+beta^mk*dk;
    new_ln=feval(fun, W, X, Y, x0); new_norm=norm(feval(gfun, W, X, Y, x0));
    if new_ln>old_ln || new_norm>old_norm
        x = x0-beta^mk*dk;
        return;
    end
    
    
    k=k+1;
end
x=x0;
val=feval(fun,W, X, Y, x0);
end