classdef Curved < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        thickness
        length
        height
        width
        radius
        sweep_angle
        start_angle
        end_angle
        youngs_mod
        shear_mod
        shape_factor
        rc_ratio
        ec_ratio
        couple

        ba % boundary angle
        bV % boundary V def (v_def(1))
        bH % boundary H def (h_def(1))

        U % strain energy
        dh % deflection parallel to cross-section
        dv % deflection through cross-sectopm
        dth % angular deflection
        Mth % sum of moments
        Vth % sum of transverse shear stresses
        Nth % sum of axial stresses
    end
    
    methods
        function obj = Curved(const, vals)
            %Curved class to model deflection of large (R>10h) curved
            %beams using Castigliano's Theorem.
            
            syms R V H Pv Ph Mo Mp th g
            syms E I G e A F
            syms l1 l2
            syms Mth(V, H, Pv, Ph, Mo, Mp, R, g, th)
            syms Vth(V, H, Pv, Ph, th)
            syms Nth(V, H, Pv, Ph, th)
            
            % Assign values to properties
            obj.youngs_mod = const('youngs_mod');
            obj.shear_mod = const('shear_mod');
            obj.shape_factor = const('shape_factor');
            
            obj.thickness = vals('thickness');

            obj.width = vals('width');
            obj.radius = vals('radius');
            obj.sweep_angle = vals('sweep_angle');
            obj.couple = vals('couple');
            obj.start_angle = vals('start_angle');
            obj.end_angle = vals('end_angle');

            obj.length = obj.radius*(cos(obj.start_angle)-cos(obj.end_angle));
            obj.height = obj.radius*(sin(obj.start_angle)-sin(obj.end_angle));
            
            obj.ba = vals('boundary_angle');
            obj.bV = vals('boundary_defV');
            obj.bH = vals('boundary_defH');

            % from 'Roark's Formulas for Stress and Strain', 9th Edition
            % for rectangular cross-sections
            obj.rc_ratio = [1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0];
            obj.ec_ratio = [0.366, 0.284, 0.236, 0.204, 0.18, 0.115, 0.085, 0.056, 0.042, 0.033];

%             % Sum of moments at angular position
%             Mth(V, H, Pv, Ph, Mo, Mv, R, g, th) = V*R*sin(th) + H*R*(1-cos(th)) + Mo + Pv*R*(2*cos((th+g)/2)*sin((th-g)/2)) - Ph*R*(2*sin((g+th)/2)*sin((g-th)/2)) + Mv;
%             Vth(V, H, Pv, Ph, th) = V*cos(th) + H*sin(th) + Pv*cos(th) + Ph*sin(th);
%             Nth(V, H, Pv, Ph, th) = -H*cos(th) - Ph*cos(th) + V*sin(th) + Pv*sin(th);

            % Sum of moments at angular position
            Mth(V, H, Pv, Ph, Mo, Mp, R, g, th) = H*R*sin(th) + V*R*(1-cos(th)) + Mo + Ph*R*(2*cos((th+g)/2)*sin((th-g)/2)) - Pv*R*(2*sin((g+th)/2)*sin((g-th)/2)) + Mp;
            
            Vth(V, H, Pv, Ph, th) = H*cos(th) + V*sin(th) + Ph*cos(th) + Pv*sin(th);
            Nth(V, H, Pv, Ph, th) = -V*cos(th) - Pv*cos(th) + H*sin(th) + Ph*sin(th);
            
            obj.Mth = Mth(V, H, Pv, Ph, Mo, Mp, R, g, th);
            obj.Vth = Vth(V, H, Pv, Ph, th);
            obj.Nth = Nth(V, H, Pv, Ph, th);

            if obj.radius/(obj.thickness/2) > 10.0
                dU = (R*Mth(V, H, Pv, Ph, Mo, Mp, R, g, th)^2)/(2*E*I);
            else
                disp('chonk')
               dU = (Mth(V, H, Pv, Ph, Mo, Mp, R, g, th)^2)/(2*A*E*e) ... % no *R in Rdth because it cancels with denominator
                + (F*(Vth(V, H, Pv, Ph, th)^2)*R)/(2*A*G) ...
                + ((Nth(V, H, Pv, Ph, th)^2)*R)/(2*A*E) ...
                - (Mth(V, H, Pv, Ph, Mo, Mp, R, g, th)*Nth(V, H, Pv, Ph, th))/(A*E);
            end

            obj.U = int(dU, th, [l1, l2]);
            obj.dv = int(diff(dU, Pv), th, [l1, l2]);
            obj.dh = int(diff(dU, Ph), th, [l1, l2]);
            obj.dth = int(diff(dU, Mp), th, [l1, l2]);
        end
        
        function [sH, sV, Hta, Vta, M, ang, U] = def(obj, forces, out_type)
            % break f into global V and H. Use to find couple and pass into object
            % f = matrix [V; H] Transform V and H locally in object
% 
%             fta_clockwise = [cos(obj.ba) sin(obj.ba); ...
%                 -sin(obj.ba) cos(obj.ba)]*forces('f'); 
            fta_cclockwise = [cos(0) -sin(0); ...
                sin(0) cos(0)]*forces('f');

%             fta_cclockwise = [cos(obj.ba) -sin(obj.ba); ...
%                 sin(obj.ba) cos(obj.ba)]*forces('f');


            % obj.ba: boundary angular deflection
            % net force already distributed to vertical and horizontal in
            % script - here we distribute the components to align
            % with the axes of the beam, deflected from bending in the
            % previous section. Vectors rotate counter clockwise with
            % positive angular deflection

            V = fta_cclockwise(1);
            H = fta_cclockwise(2);
            
            % flip forces if angles change counter clockwise
            if (obj.start_angle > obj.end_angle) % || (obj.start_angle < 0 || obj.end_angle < 0)
                H = -1*H;
                V = -1*V;
            end

            R = obj.radius;
            E = obj.youngs_mod;
            G = obj.shear_mod;
            F = obj.shape_factor;
            I = (obj.width*obj.thickness^3)/12;
            A = obj.width*obj.thickness;
            e = (obj.thickness/2)*interp1(obj.rc_ratio, obj.ec_ratio, obj.radius/(obj.thickness/2));
            
            % as these are probing forces at the point of interest
            Pv = 0;
            Ph = 0;
            Mp = 0;
            
            % subtract moment from 'ghost' first section of curve in case
            % beam does not start at 0 deg
            Mo = 0;
            g = obj.start_angle;
            th = obj.start_angle;

            Mo = obj.couple - vpa(subs(obj.Mth));

            % flip forces if angles change counter clockwise
%             if (obj.start_angle > obj.end_angle) % || (obj.start_angle < 0 || obj.end_angle < 0)
%                 Mo = -1*Mo;
%             end

            if(matches('sum', out_type))
                l1 = obj.start_angle;
                l2 = obj.end_angle;

                th = obj.start_angle; % for finding Mth
                g = obj.start_angle;
                
                % undeformed position of loaded end
                sH = R*(cos(l1))-R*(cos(l2));
                sV = R*(sin(l1))-R*(sin(l2));
                
                % clockwise roatation matrix
                rm = [cos(obj.ba) -sin(obj.ba) 0; sin(obj.ba) cos(obj.ba) 0; 0 0 1]; % clockwise
                
                % shift such that beam is roated about fixed end
                tr0 = [1 0 -0; 0 1 -0; 0 0 1];
%                 tr1 = [1 0 R*(cos(l2)); 0 1 R*(sin(l2)); 0 0 1];
                tr1 = [1 0 0; 0 1 0; 0 0 1];

                ta = [sH+vpa(subs(obj.dh)); sV+vpa(subs(obj.dv)); ones(size(sH))];
                tra = tr1*rm*tr0*ta;
                % for 'sum', return abs deformation instead of relative to static beam
                Hta = tra(1, 1:end)-sH;
                Vta = tra(2, 1:end)-sV;

                ang = vpa(subs(obj.dth));

%                 % flip back for energy, moment, angle calculations
                if (obj.start_angle > obj.end_angle) % || (obj.start_angle < 0 || obj.end_angle < 0)
%                     H = -1*H;
%                     V = -1*V;
                    M = -1*vpa(subs(obj.Mth));
                else
                    M = vpa(subs(obj.Mth));
                end
                U = abs(vpa(subs(obj.U)));

            elseif(matches('linspace', out_type))
                l1 = linspace(obj.start_angle, obj.end_angle, abs(round(R*(obj.end_angle-obj.start_angle)*1000)));
                l2 = ones(1,length(l1))*obj.end_angle;
                th = l1; % for moment obj.calc (lim from l1, l2 otherwise)
                g = l1;

                sH = R*(cos(l1));
                sV = R*(sin(l1));

                rm = [cos(obj.ba) -sin(obj.ba) 0; sin(obj.ba) cos(obj.ba) 0; 0 0 1]; % cclockwise

%                 tm1 = [1 0 -sV(end); 0 1 0; 0 0 1]; % translate to coord system origin
%                 tm01 = [1 0 sV(end); 0 1 0; 0 0 1];
% 
%                 tm2 = [1 0 -sH(end); 0 1 0; 0 0 1];
%                 tm02 = [1 0 sH(end); 0 1 0; 0 0 1];
% 
%                 Hta = tm01*rm*tm1*[sV; vpa(subs(obj.dh)); ones(size(sV))];
%                 Vta = tm02*rm*tm2*[sH; vpa(subs(obj.dv)); ones(size(sH))];
% 
%                 Hta = sH+Hta(2, 1:end);
%                 Vta = sV+Vta(2, 1:end);

                tr0 = [1 0 -sH(end); 0 1 -sV(end); 0 0 1];
                tr1 = [1 0 sH(end); 0 1 sV(end); 0 0 1];

                ta = [sH+vpa(subs(obj.dh)); sV+vpa(subs(obj.dv)); ones(size(sH))];
                tra = tr1*rm*tr0*ta;
                Hta = tra(1, 1:end);
                Vta = tra(2, 1:end);

                ang = vpa(subs(obj.dth)); % may need to rotate angle

%                 % flip back for energy, moment, angle calculations
                if (obj.start_angle > obj.end_angle) % || (obj.start_angle < 0 || obj.end_angle < 0)
%                     H = -1*H;
%                     V = -1*V;
                    M = -1*vpa(subs(obj.Mth));
                else
                    M = vpa(subs(obj.Mth));
                end
                U = abs(vpa(subs(obj.U)));
%                 M = vpa(subs(obj.Mth)); % may need to adjust moment based on change in length
                
                sH = fliplr(sH);
                sV = fliplr(sV);
                Hta = fliplr(Hta);
                Vta = fliplr(Vta);
                M = fliplr(M);
                ang = fliplr(ang);
                U = fliplr(U);
                
            else
                throw(MException('MATLAB:InvalidInputType','choose either <sum> or <linspace>'))
            end
        end
    end
end
