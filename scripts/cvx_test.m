m = 10; p = 4;
A = randn(p,m); b = randn(p,1);
I = {[1:5],[6:10]};
cvx_begin
    variable x(m)
    minimize( 0 )
    subject to
        A * x == b
        for i = 1:length(I)
            norm( x(I{i}) ) <= 1
        end
cvx_end