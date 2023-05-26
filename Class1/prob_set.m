function [c,r,l,p,q,gama] = prob_set(m,n,prob)
switch prob
    case 0  % half l1-norm = 0.5*norm(r-l,1)
        if m~=n
            error('half l1-norm problem requires m == n');
        end
        C = ones(n) - eye(n);c = C(:);
        r = rand(n,1);q = ones(n,1);
        l = rand(m,1);p = ones(m,1);
        Gama = + inf(m,n);gama = Gama(:);
        if r'*q > l'*p
            l = l + (r'*q - l'*p)/sum(p);
        else
            r = r + (l'*p - r'*q)/sum(q);
        end
    case 1  % Assignment problem
        if m~=n
            error('Assignment problem requires m == n');
        end
        C = rand(m,n); c = C(:);
        r = ones(n,1);q = ones(n,1);
        l = ones(m,1);p = ones(m,1);
        Gama = + inf(m,n);gama = Gama(:);
    case 2  % Optimal transport
        C = rand(m,n); %example 1
%         C = ones(m,n)-eye(m,n); %example 2
%         cx = rand(m,1); cy = rand(n,1); C = (cx*ones(1,n)-ones(m,1)*cy').^2; %example 3 
        c = C(:);
        r = rand(n,1);q = ones(n,1);
        l = rand(m,1);p = ones(m,1);
        Gama = + inf(m,n);gama = Gama(:);
        if r'*q > l'*p
            l = l + (r'*q - l'*p)/sum(p);
        else
            r = r + (l'*p - r'*q)/sum(q);
        end
    case 3  % Capacity constrained optimal transport
%         C = rand(m,n); %example 1
%         C = ones(m,n)-eye(m,n); %example 2
        cx = rand(m,1); cy = rand(n,1); C = (cx*ones(1,n)-ones(m,1)*cy').^2; %example 3
        c = C(:);
        r = rand(n,1);q = ones(n,1);
        l = rand(m,1);p = ones(m,1);
        if r'*q > l'*p
            l = l + (r'*q - l'*p)/sum(p);
        else
            r = r + (l'*p - r'*q)/sum(q);
        end
        Gama = l*r'/sum(l)+10*rand(m,n);gama = Gama(:);
end
end