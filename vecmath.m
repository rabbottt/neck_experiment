classdef vecmath
    
    methods (Static)
        
        % Transform data to new scale based on [min max] 
        function new_data = changeScale(data, iR, oR)
            % Inputs:       iR: (nx2) input Range n rows of [min max] 
            %               oR: (nx2) output Range n rows of [min max]
            %               data: (nx1) vector of data in input range (corresponding to [min max] in nth row)
            % Outputs:      new_data: (nx1) vector of data in output range
            % 
                       
            new_data = (oR(:,2)-oR(:,1))./(iR(:,2)-iR(:,1)).*(data - iR(:,1)) + oR(:,1);    
        end
        
        function maxn = maxnorm(ranges)
            % ranges: (n x [min max]) matrix 
            % maxn: (scaler) maximum norm possible from range input
            maxn = norm(max(abs(ranges),[],2));
        end
        
        function [uv, mag] = unitvec(v)
            % v: any vector
            % uv: unit vector of v
           
            % compute magnitude
            mag = norm(v);
            
            % compute unit vector
            if mag ~= 0 
                uv = v/mag; 
            else % mag == 0, avoid division by 0
                uv = 0*v;
            end
        end
        
        function th = aBetween2Vecs(v1,v2)
            % Angle between any 2 vectors of same length
            % v1, v2: vectors of same dimensions
            % th: angle (radians) between v1 and v2
            %     if v1 or v2 are zero vectors, th = NaN
            
            % magnitudes
            v1mag = norm(v1);
            v2mag = norm(v2);
            
            % compute angle between v1 and v2
            if any([v1mag == 0, v2mag == 0]) % avoid division by 0
                th = NaN; % not valid
            else % valid           
                th = acos(dot(v1,v2)/(v1mag*v2mag));
            end

        end
        
        % (3x3) tranformation matrix to:
        %   reflect about a plane
        function R = ReflectAboutPlane(plane)
            % Inputs:       plane: 'xy', 'yz', or 'xz'
            % Outputs:      R : 3x3 matrix to reflect about plane 
            
            switch plane
                case 'xy'
                    dvec = [1 1 -1];
                case 'yz'
                    dvec = [-1 1 1];
                case 'xz'
                    dvec = [1 -1 1];
                otherwise 
                    dvec = [1 1 1];
                    warning('options for plane are xy yz and xz')
            end

            R = diag(dvec);
            
        end
        
        
        % (3x3) transformation matrix to:
        % rotate around arbitrary axis u (1x3) by th (scaler) in radians
        function R = RotateAxis(u, th)
            % u: unit vector, arbitrary axis to rotate around
            % th: angle in radians to rotate about u
            % R : 3x3 rot mat
            
            R = eye(3);
            nu = norm(u);
            eps = 0.01; % sensitivity for detecting rotation direction 
            if nu > eps % if not a zero vector
                % make sure it is a unit vector
                u = u/nu;
                % trig
                s = sin(th);
                c = cos(th);
                % Rodriguez formula
                for ii=1:3
                    v = R(:,ii);
                    R(:,ii) = v*c + cross(u,v)*s + u*dot(u,v)*(1-c);                      
                end                
            end          
        end
        
        function R = RotateVec2Vec(vec1, vec2)
            % 3D rotation matrix to align vec1 to vec2            
            
            % make sure they are unit vectors
            v1 = vec1/norm(vec1);
            v2 = vec2/norm(vec2);            
            
            if v1 == v2 
                % identity
                R = eye(3);
            else       
                % 
                v = cross(v1,v2);
                
                % skew-symmetric matrix 
                ssc = [ 0       -v(3)    v(2); ...
                        v(3)     0      -v(1); ...
                       -v(2)     v(1)    0  ];
                % compute rotation matrix 
                R = eye(3) + ssc + ssc^2*(1-dot(v1,v2))/(norm(v)^2);
            end
            
        end
        
        function Rx = RotX(th)
            Rx = [  1   0        0;...
                    0   cos(th) -sin(th);...
                    0   sin(th)  cos(th)];
        end

        function Ry = RotY(th)
            Ry = [  cos(th)     0   sin(th)   ;...
                    0           1   0         ;...
                    -sin(th)    0   cos(th)   ];        
        end

        function Rz = RotZ(th)
            Rz = [  cos(th)     -sin(th)    0   ;...
                    sin(th)      cos(th)    0   ;...
                    0            0          1   ];
        end
    
        % homogenous (4x4) transformation matrix to: 
        % rotate by R (3x3) and translate by p (1x3) or (3x1)
        function Tmat = TransformMat(R,p)
            Tmat = [R [p(1) p(2) p(3)]'; [0 0 0 1]];
        end
    
        % Adjoint matrix to transform wrench coordinate frames
        function adT = adjointMatrix(R,p)
            % R: (3x3) rotation matrix 
            % p: (3x1) or (1x3) translation vector
            % adT: (6x6) adjoint matrix 
            % Wb = adT' * Wa % to change coordinates from a to b frame
            % W must be in form [F M]' NOT [M F]
            
            % skew-symmetric matrix 
            phat = [    0    -p(3)  p(2); ...
                        p(3)  0    -p(1); ...
                       -p(2)  p(1)  0  ];
        
            adT = [ R           phat*R   ; ...
                    zeros(3,3)  R       ];
        end
        
        % homogenous (4x4) transformation matrix to:
        % scale by s (scalar) uniformly about origin
        function Tmat = scaleMat(s)            
            Tmat = [s*eye(3) zeros(3,1); [0 0 0 1]];            
        end
        
        function R = quat2R(qvec)
        % Transform Quaternion into 3x3 Rotation Matrix 
        % IN: 
        %     Qvec - quaternion describing rotation
        % 
        % OUT:
        %     R - rotation matrix 

            % Break down Quaternion into components
            w = qvec(1);
            x = qvec(2);
            y = qvec(3);
            z = qvec(4);
            
            % Rotation Matrix Components
            Rxx = 1 - 2*(y^2 + z^2);
            Rxy = 2*(x*y - z*w);
            Rxz = 2*(x*z + y*w);
            Ryx = 2*(x*y + z*w);
            Ryy = 1 - 2*(x^2 + z^2);
            Ryz = 2*(y*z - x*w );
            Rzx = 2*(x*z - y*w );
            Rzy = 2*(y*z + x*w );
            Rzz = 1 - 2 *(x^2 + y^2);
            
            % Build Rotation Matrix
            R = [   Rxx,    Rxy,    Rxz;
                    Ryx,    Ryy,    Ryz;
                    Rzx,    Rzy,    Rzz];            
            
        end
    
        function pts = hom(pts) 
            % adds additional ROW of ones
            npts = size(pts,2);
            pts = [pts; ones(1,npts)];
        end
        
        function pts = unhom(pts)
            % removes last ROW
            pts = pts(1:3,:);
        end
        
        function [v,f] = createRect(x0,y0,xm,ym)
            % vertices for rectangle 
            v = [  x0      y0; ... % 1 % lower left
                   xm      y0; ... % 2 % lower right
                   xm      ym; ... % 3 % upper right
                   x0      ym];    % 4 % upper left
            
            f = [ 1 2 3 4 ];       % face
        end
        
        % return bounded value clipped between lb and ub
        function y = clip(x,lb,ub)
            y = min( max(x, lb), ub);
        end  
        
        
    end
end