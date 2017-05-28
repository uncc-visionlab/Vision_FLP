classdef FittingLibrary < handle
    methods (Static)
        function [coeffs, error] = fitPlane(pts, subtract_mean)
            % Fits an implicit plane equation to a list of points
            % Each point should be on a different row.
            % there should be at least as many points as there are dimensions.
            if(exist('subtract_mean')==0)
                subtract_mean=1;
            end
            if (size(pts,1) < size(pts,2))
                pts=pts';
            end
            if (subtract_mean==1)
                meanVal=mean(pts);
                zpts=pts-ones(size(pts,1),1)*meanVal;
                [evec,eval]=eig(zpts'*zpts);
                coeffs = evec(:,1); % vector associated with smallest eigenvalue
                d = -meanVal*coeffs;
                normf = norm(coeffs(1:3));
                error_sq = eval(1,1);
                error = sqrt(error_sq)/(normf*size(pts,1));
                coeffs = [coeffs' d]./normf;
            else
                npts=[pts ones(size(pts,1),1)];
                npts'*npts
                [evec,eval]=eig(npts'*npts);
                coeffs = evec(:,1); % vector associated with smallest eigenvalue
                normf = norm(coeffs(1:3));
                coeffs = coeffs'./normf;
                error_sq = eval(1,1);
                error = sqrt(error_sq)/(normf*size(pts,1));
            end
        end
        
        function [coeffs, error] = fastFitPlane(theta_x, theta_y, z)
            % Fits an implicit plane equation to a list of points
            % Each point should be on a different row.
            % there should be at least as many points as there are dimensions.
            
            M=[theta_x theta_y ones(size(theta_x,1),1) 1./z];
            %meanVal=mean(M);
            [evec,eval]=eig(M'*M);
            coeffs = evec(:,1);  % vector associated with smallest eigenvalue
            normf = norm(coeffs(1:3));
            coeffs = coeffs'./normf;
            error_sq = eval(1,1);
            error = sqrt(error_sq)/(normf*size(theta_x,1));
            
            if (1==0)
                M1=[theta_x theta_y ones(size(theta_x,1),1)];
                %meanVal=mean(M);
                coeffs2 = inv(M1'*M1)*M1'*(1./z);
                coeffs2 = [coeffs2; -1];
                normf = norm(coeffs2(1:3));
                coeffs2 = coeffs2'./normf
                coeffs=coeffs2;
            end
            %error_sq = eval(1,1);
            %error = sqrt(error_sq)/(normf*size(theta_x,1));
            %errors=M*coeffs';
            %error=sqrt(errors'*errors)/size(theta_x,1);
            %error-sqrt(error_sq)/size(theta_x,1)
        end
        
        function [coeffs, error] = fitQuadric(pts)
            if (size(pts,1) < size(pts,2))
                pts=pts';
            end
            N = size(pts,1);
            
            ptcloud_mean = mean(pts);
            pts=pts-ones(size(pts,1),1)*ptcloud_mean;
            aabb=[min(pts); max(pts)];
            
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            M=[x.^2 x.*y x.*z y.^2 y.*z z.^2 x y z ones(N,1)];
            Gx=0.001*[2.*x y z zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1)];
            Gy=0.001*[zeros(N,1) x zeros(N,1) 2.*y z zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1)];
            Gz=0.001*[zeros(N,1) zeros(N,1) x zeros(N,1) y 2.*z zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1)];
            A=[M;Gx;Gy;Gz];
            b=[zeros(N,1); ones(N,1); ones(N,1); ones(N,1)];
            coeffs2=pinv(A)*b;
            
            S=M'*M;
            %S(7:10,7:10)
            [evec,eval]=eig(S(7:10,7:10));
            plane_coeffs=evec(:,1);
            normf=norm(plane_coeffs(1:3));
            plane_coeffs = plane_coeffs'/normf;
            if (plane_coeffs(3) < 0)
                plane_coeffs = -plane_coeffs;
            end
            %plane_coeffs
            plane_error=sqrt(eval(1,1))/(normf*N);
            [evec,eval]=eig(S);
            quad_coeffs=evec(:,1);
            quad_error=sqrt(eval(1,1))/(normf*N);
            quad_coeffs = coeffs2;
            % df_dx = 2.*quad_coeffs(1).*x + quad_coeffs(2).*y + quad_coeffs(3).*z + quad_coeffs(7).*ones(N,1);
            % df_dy = quad_coeffs(2).*x + 2.*quad_coeffs(4).*y + quad_coeffs(5).*z + quad_coeffs(8).*ones(N,1);
            % df_dz = quad_coeffs(3).*x + quad_coeffs(5).*y + 2.*quad_coeffs(6).*z + quad_coeffs(9).*ones(N,1);
            % gradMagnitude=sqrt(df_dx.^2+df_dy.^2+df_dz.^2);
            % df_dx = df_dx./gradMagnitude;
            % df_dy = df_dy./gradMagnitude;
            % df_dz = df_dz./gradMagnitude;
            % normVecMat=[df_dx df_dy df_dz];
            %
            % normVecMat(normVecMat(:,3)<0,:)=-normVecMat(normVecMat(:,3)<0,:) % delete rows where depth is Z=0
            
            cc = quad_coeffs;
            Q=0.5*[0 cc(2) cc(3);
                cc(2) 0 cc(5);
                cc(3) cc(5) 0];
            Q(1,1)=cc(1);
            Q(2,2)=cc(4);
            Q(3,3)=cc(6);
            P(1)=cc(7);
            P(2)=cc(8);
            P(3)=cc(9);
            [R1,V,R2]=svd(Q)
            
            qqq=R2'*Q*R2
            ppp2=P*R2
            translation_val = ppp2./(2*diag(qqq)')
            constant_mod = (1/4)*((ppp2.^2)./diag(qqq)');
            
            coeffs_new=[qqq(1,1) 0 0 qqq(2,2) 0 qqq(3,3) 0 0 0 cc(10)-sum(constant_mod)];
            
            cc = coeffs_new;
            xyz_new = (R2'*pts')';
            xyz_new = xyz_new+ones(size(xyz_new,1),1)*translation_val;
            aabb_new = (R2'*aabb')'+ones(size(1,1),1)*translation_val;
            aabb_check=[min(xyz_new); max(xyz_new)];
            aabb_new = aabb_check;
            % cc(6)=0;
            % cc(10)=0;
            x=xyz_new(:,1);
            y=xyz_new(:,2);
            z=xyz_new(:,3);
            df_dx = 2.*cc(1).*x + cc(2).*y + cc(3).*z + cc(7).*ones(N,1);
            df_dy = cc(2).*x + 2.*cc(4).*y + cc(5).*z + cc(8).*ones(N,1);
            df_dz = cc(3).*x + cc(5).*y + 2.*cc(6).*z + cc(9).*ones(N,1);
            %gradVals=[df_dx df_dy df_dz];
            normMags=sqrt(df_dx.^2+df_dy.^2+df_dz.^2);
            df_dx = df_dx./normMags;
            df_dy = df_dy./normMags;
            df_dz = df_dz./normMags;
            normVals=[df_dx df_dy df_dz];
            figure(3), hold off, plot3(xyz_new(:,1), ...
                xyz_new(:,2),xyz_new(:,3),'b.','MarkerSize',12);
            
            figure(3), hold on, quiver3(xyz_new(:,1),xyz_new(:,2),xyz_new(:,3), ...
                normVals(:,1),normVals(:,2),normVals(:,3),'r');
            
            %aabb_new(1,:)=max(xyz_new);
            %aabb_new(2,:)=min(xyz_new);
            npts=100;
            [y,x,z] = ndgrid(linspace(aabb_new(1,2),aabb_new(2,2),npts), ...
                linspace(aabb_new(1,1),aabb_new(2,1),npts), ...
                linspace(aabb_new(1,3),aabb_new(2,3),npts));
            f = cc(1).*x.^2 + cc(2).*x.*y + cc(3).*x.*z  + cc(4).*y.^2 + ...
                cc(5).*y.*z + cc(6).*z.^2 + cc(7).*x + cc(8).*y + cc(9).*z + cc(10);
            %cc(6)=0;
            %cc(7:10)=0;
            cc(6)=0;
            cc(10)=0;
            f2 = cc(1).*x.^2 + cc(2).*x.*y + cc(3).*x.*z  + cc(4).*y.^2 + ...
                cc(5).*y.*z + cc(6).*z.^2 + cc(7).*x + cc(8).*y + cc(9).*z + cc(10);
            cc
            
            figure(3), hold on, isosurface(x,y,z,f,.01);
            % figure(3), hold on, isosurface(x,y,z,f2,.00001);
            view(3);
            camlight
            axis equal
            
            line_segment = [0 0 aabb_new(1,3); 0 0 aabb_new(2,3)];
            figure(3), plot3(line_segment(:,1),line_segment(:,2),line_segment(:,3),'g');
            % line_segment = line_segment-ones(2,1)*translation_val;
            % figure(2), plot3(line_segment(:,1),line_segment(:,2),line_segment(:,3),'g');
            % line_segment = (R2*line_segment')'
            % figure(1), plot3(line_segment(:,1),line_segment(:,2),line_segment(:,3),'g');
            
            
            %avg_angle=(180/pi)*sum(abs(acos(normVecMat*plane_coeffs(1:3)')))/N;
            %meanVec=mean(normVecMat)
            %cov(normVecMat)
            coeffs=plane_coeffs;
            %error=avg_angle;
            error=plane_error;
        end
    end
end