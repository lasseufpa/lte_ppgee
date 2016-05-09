function app = bcjr2(trellis,priori, channel)
priori
channel
k = log2(trellis.numInputSymbols);
n = log2(trellis.numOutputSymbols);

numIn = trellis.numInputSymbols;
numSt = trellis.numStates;
next = trellis.nextStates;
output = trellis.outputs;

prev = zeros(numSt, numIn);
for s = 1:numSt
   for i=1:numIn 
       nextSt = next(s,i);
       prev(nextSt+1, i) = s-1;
   end
end

l = length(channel)/n;

A = zeros(l+1, numSt);
B = zeros(l+1, numSt);
L = zeros(l+1, numSt);
G = zeros(l, numSt, numIn);
M = zeros(l, numSt, numIn);
app = zeros(l*k,1);
for t=1:l
    for s=1:numSt
        for i=1:numIn
            G(t,s,i) = 1;
            for j=1:k
                p = priori((t-1)*k+j);
                u = bitget(i-1,k-j+1);
                if(u==1)
                    p = 1-p;
                end
                G(t,s,i) = G(t,s,i) * p;
            end
            o = output(s,i);
            for j=1:n
                p = channel((t-1)*n+j);
                u = bitget(o,n-j+1);
                if(u==1)
                    p = 1-p;
                end
                G(t,s,i) = G(t,s,i) * p;
            end
        end
    end
end
G
A(1,1) = 1;
for t=2:l+1
    for s = 1:numSt
        for i = 1:numIn
            prevSt = prev(s,i)+1;
            A(t,s) = A(t,s) + A(t-1, prevSt) * G(t-1,prevSt, i); 
        end
    end
end

B(l+1,:) = 1/numSt;
for t=l:-1:1
    for s = 1:numSt
        for i=1:numIn
            nextSt = next(s,i)+1;
            B(t,s) = B(t,s) + B(t+1,nextSt) * G(t, s, i);
        end
    end
end

for t=1:l+1
    for s=1:numSt
        for i=1:numIn
            L(t,s) = A(t,s)*B(t,s);
        end
    end
end
norm = sum(L,2);
for t=1:l+1
    L(t,:) = L(t,:) ./ norm(t);
end
L

for t=1:l
    for s=1:numSt
        for i=1:numIn
            nextSt = next(s,i)+1;
            M(t,s,i) = A(t,s) * G(t,s,i) * B(t+1,nextSt);
        end
    end
end
norm = sum(sum(M,2),3);
for t=1:l
    M(t,:,:) = M(t,:,:)./norm(t);
end
M

for t=1:l
    for j=1:k
        p0 = 0;
        for s=1:numSt
            for i=1:numIn
                u = bitget(i-1,k-j+1);
                if(u==0)
                    p0 = p0 + M(t,s,i);
                end
            end
        end
        app((t-1)*k+j)=p0;
    end
end

end