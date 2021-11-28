function J=Yang_Retinex(I)  
%%%%%%%%%%Define parameters
        Sx=0.0023;
        Lamdar=620;
        Lamdag=540;
        Lamdab=450;
        krnlsz=5;
        r = krnlsz*6;
        eps = 10^-6;


        %%%%%%%%%%%%% Compute number of channel of input image
        [m,n,z]=size(I);
        R=I(:,:,1);
        G=I(:,:,2);
        B=I(:,:,3);
%         I_gray=rgb2gray(I);
        %%If the dark channel hypothesis is met, estimate tr
        %%The image block size krnlsz is defined and DC is obtained
                tic
        DC=Dark_method(I,15);
            DC=Dark_Channel_Red(I);
%             DC2=Dark_Channel_GB(I);

        %%%%The second method is to find bwq
        sigma1=10;%double(uint8(h*0.03));%The small size is 3% of the image size.
        F1=fspecial('gaussian',[15,15],sigma1);%Small size template.
        Bwqr=imfilter(R,F1,'replicate','conv');%The three Gaussian templates are convolved with the R channel.
        Bwqg=imfilter(G,F1,'replicate','conv');%The three Gaussian templates are convolved with the G channel.
        Bwqb=imfilter(B,F1,'replicate','conv');%The three Gaussian templates are convolved with the B channel.
        %
        % [Bwqr,Bwqg, Bwqb] = backlight(I,DC);
        % [Bwqr,Bwqg,Bwqb]=Bwq_Cosmin(I,DC);
        %         Bwqr(Bwqr>1)=1;
        %         Bwqg(Bwqg>1)=1;
        %         Bwqb(Bwqb>1)=1;
        % tr=1-DC;
        % DC_Red=Dark_Channel_Red(I);
        %  t_Red=1-DC_Red;
        %  tr=t_Red;
        %     tr=1-retinex_mccann99(DC, 2);
        [tr]=1-variance_retinexfu(DC);
        %                 tr=Ngretinex(DC);
        %                 tr=1-Kimmelretinex(DC.*255);
        %                               tr =retinex_central_gray(DC);
        tr=max(tr,0.1);
        %         tr(tr>1)=1;
%         tr= (guidedfilter(R, tr, r, eps));
        t_R=tr;

        %%The first way to estimate tg,tb
        kr=exp(-Sx*(Lamdag-440))/exp(-Sx*(Lamdar-440));
        kb=exp(-Sx*(Lamdab-440))/exp(-Sx*(Lamdar-440));
        t_G=real(t_R.^kr);
        t_B=real(t_R.^kb);

        %         t_G=(guidedfilter(G, t_G, r, eps));
        %         t_B=(guidedfilter(B, t_B, r, eps));
        t_G=max(t_G,0.1);
        t_B=max(t_B,0.1);

        %%The second way to estimate TG, TB
        %%5.find Jr',Jg',K=tg/tr;
        %     trp=1./tr;
        %     [dx,dy]=gradient(R);
        %     Jr=sqrt(dx.^2+dy.^2);
        %     [dx,dy]=gradient(G);
        %     Jg=sqrt(dx.^2+dy.^2);
        %     [dx,dy]=gradient(B);
        %     Jb=sqrt(dx.^2+dy.^2);
        %     %     Jrg=minfilt2(Jr./(Jg+eps), [krnlsz,krnlsz]);
        %     %   Jrg(m,n)=0;
        %     Jrg=Jr./(Jg+eps);
        %     %     Jrb=minfilt2(Jr./(Jb+eps), [krnlsz,krnlsz]);
        %     %     Jrb(m,n)=0;
        %     Jrb=Jr./(Jb+eps);
        %     K1=trp./(Jrg.*(trp-1)+1);
        %     tg=K1.*tr;
        %     K2=trp./(Jrb.*(trp-1)+1);
        %     tb=K2.*tr;
        %     tg= (guidedfilter(G, tg, r, eps));
        %     tb= (guidedfilter(B, tb, r, eps));
        %     %     tg(tg>1)=1;
        %     %     tb(tb>1)=1;

        %%Various methods for estimating Bwq
        %%5.find tr',tg',K=tg/tr;
        %     [dx,dy]=gradient(tr);
        %     trp=1./(sqrt(dx.^2+dy.^2)+eps);
        %     [dx,dy]=gradient(tg);
        %     tgp=1./(sqrt(dx.^2+dy.^2)+eps);
        %     [dx,dy]=gradient(tb);
        %     tbp=1./(sqrt(dx.^2+dy.^2)+eps);
        %     Bwqr=Jr.*trp;
        %     Bwqg=Jg.*tgp;
        %     Bwqb=Jb.*tbp;

        %     t_B=max(t_B,0.1);
        %     t_G=max(t_G,0.1);

        J = zeros(m,n,3);
        J(:,:,1) = (I(:,:,1)-Bwqr)./t_R+Bwqr;
        %         J(:,:,1) = (J(:,:,1) - min(min(J(:,:,1))))/(max(max(J(:,:,1))) -  min(min(J(:,:,1))));
        J(:,:,2) = (I(:,:,2)-Bwqg)./t_G+Bwqg;
        %         J(:,:,2) = (J(:,:,2) - min(min(J(:,:,2))))/(max(max(J(:,:,2))) -  min(min(J(:,:,2))));
        J(:,:,3) = (I(:,:,3)-Bwqb)./t_B+Bwqb;
        %         J(:,:,3) = (J(:,:,3) - min(min(J(:,:,3))))/(max(max(J(:,:,3))) -  min(min(J(:,:,3))));
        %         J=gammajiaozheng(J.*255,5);


        for i=1:z

            u=mean2(J(:,:,i));%Calculate the mean, variance, minimum and maximum value of the image.
            s=std2(J(:,:,i));
            mi=u-2.1*s;
            ma=u+2.1*s;
            J(:,:,i)=(J(:,:,i)-mi)/(ma-mi);%Direct contrast drawing.
            J(J < 0) = 0;
            J(J > 1) = 1;
        end
%         HSV=rgb2hsv(J);
%         %         for i=2:2
%         HSV(:,:,2) = imadjust( HSV(:,:,2));
%         %         u=mean2(HSV(:,:,i));%Calculate the mean, variance, minimum and maximum value of the image.
%         %         s=std2(HSV(:,:,i));
%         %         mi=u-3*s;
%         %         ma=u+3*s;
%         %         HSV(:,:,i)=(HSV(:,:,i)-mi)/(ma-mi);%Direct contrast drawing.
%         HSV(HSV< 0) = 0;
%         HSV(HSV > 1) = 1;
%         %         end
%         %
%         J=hsv2rgb(HSV);