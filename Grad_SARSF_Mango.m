function x=Grad_SARSF_Mango(fun,gfun, x0, W, X, Y, epsilon)
%功能: 梯度法求解无约束优化问题:  min f(x)
%输入: fun, gfun分别是目标函数及其梯度, x0是初始点,
%      epsilon为容许误差
%输出: k是迭代次数, x, val分别是近似最优点和最优值
maxk=1000;  %最大迭代次数
beta=0.1;  sigma=0.4;
k=0;
while(k<maxk)
    gk=feval(gfun,W, X, Y, x0);  %计算梯度
    old_ln=feval(fun, W, X, Y, x0);
    dk=-gk;  %计算搜索方向
    old_norm=norm(gk);
    if(old_norm<epsilon), break; end  %检验终止准则
    m=0;  mk=0;
    while(m<20)  %用Armijo搜索求步长
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