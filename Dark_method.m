    function DC=Dark_method(R,krnlsz,G,B)
    [m,n,z]=size(R);
    if nargin<2
        I1=R;
        R=zeros(m,n);
        R=I1(:,:,1);
        G=I1(:,:,2);
        B=I1(:,:,3);
        krnlsz=15;
    else if nargin<3
                  I1=R;
        R=zeros(m,n);
        R=I1(:,:,1);
        G=I1(:,:,2);
        B=I1(:,:,3);
        end
    end
%   eps = 10^-6;

    MR(:,:,1)=minfilt2(R, [krnlsz,krnlsz]);
    MR(:,:,2)=minfilt2(G, [krnlsz,krnlsz]);
    MR(:,:,3)=minfilt2(B, [krnlsz,krnlsz]);
    MR(m,n,1)=0;
    MR(m,n,2)=0;
    MR(m,n,3)=0;
    for y=1:m
        for x=1:n
            DC(y,x) = min(MR(y,x,:));%Minimum of four channels
        end
    end