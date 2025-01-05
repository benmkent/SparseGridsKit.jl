import PolyChaos: clenshaw_curtis

transformdomain(xw,a,b) = (a .+ (b-a) .* 0.5 .* (xw[1].+1.0) , xw[2])

doubling(n) = n == 1 ? 1 : 2^(n - 1) + 1

linear(n) = n

linearpoints(n) = 1:n, ones(n)

uniformpoints(n, a, b) = transformdomain(uniformpoints(n),a,b)
uniformpoints(n) = n==1 ? ([0.0] , [1.0]) : (range(-1, stop=1, length=n), 1/n .* ones(n))

ccpoints(n, a, b) =  transformdomain(ccpoints(n),a,b)
function ccpoints(n)
    if n==1 
        x = [0.0]
        w = [1.0]
    else
        x,w =clenshaw_curtis(n-1)
        x = 0.5*(x.-x[end:-1:1]) # Force symmetry
        w = 0.5*w
    end
    return x,w
end