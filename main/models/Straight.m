classdef Straight < handle
    %STRAIGHT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        thickness
        length % length of beam along cosine of beam angle
        height % length of beam aloing sine of beam angle
        width 
        youngs_mod % E
        shear_mod % G
        shape_factor % 6/5 for rectangular cross-section
        beam_angle % angle from 0
        section_length % abs length of beam irrespective of angle
        couple % coupling moment from preceding stitched beam section

        U % strain energy
        dh % deflection parallel to cross-section
        dv % deflection through cross-sectopm
        dth % angular deflection in plane of curvature

        ba % boundary angle (combine with angle of beam section)
        bV % boundary V def (v_def(1))
        bH % boundary H def (h_def(1))
        ca

        Ms % sum of moments
        Vs % sum of transverse shear forces
        Ns % sum of axial forces
    end
    
    methods
        function obj = Straight(const, vals)
            %STRAIGHT class to model deflection of thin straight
            %beams using obj.castigliano's Theorem.
            
            syms V H Mo E I A G F s l
            syms Ph Pv Mp lp
            syms l1 l2
            syms Ms(H, V, Mo, Ph, Pv, Mp, lp, s, l)
            syms Vs(H, V, Ph, Pv, s)
            syms Ns(H, V, Ph, Pv, s)

            % Assign input variables to properties
            obj.youngs_mod = const('youngs_mod');
            obj.shear_mod = const('shear_mod');
            obj.shape_factor = const('shape_factor');
            
            % for deflection transform to work in 2D, beam must be in 2D
%             if vals('beam_angle') ~= 0
%                 obj.beam_angle = vals('beam_angle');
%             else
%                 obj.beam_angle = 10e-10;
%             end
            obj.thickness = vals('thickness');
            obj.width = vals('width');
            obj.beam_angle = vals('beam_angle');
            obj.length = vals('length')*cos(obj.beam_angle);
            obj.height = vals('length')*sin(obj.beam_angle);
            obj.section_length = vals('length');
            obj.couple = vals('couple');
            
            % boundary conditions from any stitched beam section
            obj.ba = vals('boundary_angle');
            obj.bV = vals('boundary_defV');
            obj.bH = vals('boundary_defH');
 
            % define moment and shear expressions throughout beam
            % s: var, lo: probe length, l: total beam length
            Ms(H, V, Mo, Ph, Pv, Mp, lp, s, l) = (V*(l-s) + Pv*(lp-s))*cos(obj.beam_angle)...
                - (H*(l-s) + Ph*(lp-s))*sin(obj.beam_angle) + Mo + Mp;

            Vs(H, V, Ph, Pv, s) = (V + Pv)*cos(obj.beam_angle) - (H + Ph)*sin(obj.beam_angle);
            Ns(H, V, Ph, Pv, s) = (H + Ph)*cos(obj.beam_angle) + (V + Pv)*sin(obj.beam_angle);
            
            % account for transverse shear and normal stress
            dU = (Ms(H, V, Mo, Ph, Pv, Mp, lp, s, l)^2)/(2*E*I)...
                + (Ns(H, V, Ph, Pv, s)^2)/(2*E*A)...
                + (F*Vs(H, V, Ph, Pv, s)^2)/(2*A*G);
            
            % partial diff according to Castigliano's 2nd theorem
            obj.U = int(dU, s, [l1 l2]);
            obj.dh = int(diff(dU, Ph), s, [l1 l2]);
            obj.dv = int(diff(dU, Pv), s, [l1 l2]);
            obj.dth = int(diff(dU, Mp), s, [l1 l2]);
            
            % assign equations to properties
            obj.Ms = Ms(H, V, Mo, Ph, Pv, Mp, lp, s, l);
            obj.Vs = Vs(H, V, Ph, Pv, s);
            obj.Ns = Ns(H, V, Ph, Pv, s);

        end
        
        function [sH, sV, Hta, Vta, M, ang, U] = def(obj, forces, out_type)
%             fta_cclockwise = [cos(obj.ba) -sin(obj.ba); ...
%                 sin(obj.ba) cos(obj.ba)]*forces('f');
            fta_cclockwise = [cos(0) sin(0); ...
                -sin(0) cos(0)]*forces('f');
            V = fta_cclockwise(1);
            H = fta_cclockwise(2);
            
            Ph = 0;
            Pv = 0;
            Mp = 0;

            % substitute values for material properties and geometry
            A = obj.thickness*obj.width;
            E = obj.youngs_mod;
            G = obj.shear_mod;
            F = obj.shape_factor;
            I = (obj.width*obj.thickness^3)/12;
            Mo = obj.couple;


            if(matches('sum', out_type)) % find deflection across full beam
                l2 = obj.section_length;
                l1 = 0;

                lp = l2;
                l = obj.section_length;

                % interested in moment across structure, not at end
                % necessarily 
                s = 0; % for moment calc (lim from l1, l2 otherwise)
                M = vpa(subs(obj.Ms));

                ang = vpa(subs(obj.dth));
                U = vpa(subs(obj.U));

                sH = l2.*cos(obj.beam_angle);
                sV = l2.*sin(obj.beam_angle);

                rm = [cos(obj.ba) -sin(obj.ba) 0; sin(obj.ba) cos(obj.ba) 0; 0 0 1]; % cc
                tr0 = [1 0 -l1*cos(obj.beam_angle); 0 1 -l1*sin(obj.beam_angle); 0 0 1];
                tr1 = [1 0 l1*cos(obj.beam_angle); 0 1 l1*sin(obj.beam_angle); 0 0 1];

                ta = [sH+vpa(subs(obj.dh)); sV+vpa(subs(obj.dv)); ones(size(sH))];
                tra = tr1*rm*tr0*ta;
                Hta = tra(1, 1:end)-sH;
                Vta = tra(2, 1:end)-sV;
               

%                 lp = l1; % loaded end is x=0
%                 s = l1;
%                 
%                 M = vpa(subs(obj.Ms)); % for x = length
%                 ang = vpa(subs(obj.dth)); % may need to rotate angle
%                 U = vpa(subs(obj.U));
% 
% %                 v = vpa(subs(obj.dv));
% %                 delta = vpa(subs(obj.dh));
% %                 h = delta*cos(obj.beam_angle);
% %                 v = delta*sin(obj.beam_angle);
% 
%                 sH = l2; %*cos(obj.beam_angle); % length perpendicular to H
%                 sV = 0; % height perpendicular to V
% 
%                 rm = [cos(obj.ca) -sin(obj.ca) 0; sin(obj.ca) cos(obj.ca) 0; 0 0 1]; % clockwise
%                 tm1 = [1 0 -l2; 0 1 0; 0 0 1];
%                 tm01 = [1 0 -obj.dV+l2*cos(obj.beam_angle); 0 1 obj.bH+l2*sin(obj.beam_angle); 0 0 1];
% 
%                 Hta = tm01*rm*tm1*[sH; vpa(subs(obj.dh)); 1];
%                 Vta = [sV; 0; 1];

            elseif(matches('linspace', out_type))
                l1 = ones(1, abs(round(obj.section_length*1000)))*0; % integrating from right to left
%                 l2 = ones(1, abs(round(obj.section_length*1000)))*obj.section_length; % integrating from right to left
%                 len_l2 = size(l2);

%                 l1 = linspace(0, obj.section_length, length(l2)); % from 'full beam' to far left (l1 = lo)
                l2 = linspace(0, obj.section_length, length(l1)); % from 'full beam' to far left (l1 = lo)
%                 l1 = linspace(obj.section_length, 0, length(l2)); % from 'full beam' to far left (l1 = lo)
                lp = l2;
                l = obj.section_length;
                s = l2; % for moment obj.calc (lim from l1, l2 otherwise)

                M = vpa(subs(obj.Ms));
                ang = vpa(subs(obj.dth));
                U = vpa(subs(obj.U));

                sH = l2.*cos(obj.beam_angle);
                sV = l2.*sin(obj.beam_angle);

                % transformation matrices

%                 tm1 = [1 0 -sH(1); 0 1 0; 0 0 1]; % translate to coord system origin (no translation)
%                 tm01 = [1 0 sH(1); 0 1 0; 0 0 1];
%                 tm2 = [1 0 -sV(1); 0 1 0; 0 0 1]; % translate to coord system origin (no translation)
%                 tm02 = [1 0 sV(1); 0 1 0; 0 0 1];
%                 Hta = tm01*rm*tm1*[sV; vpa(subs(obj.dh)); ones(size(sH))];
%                 Vta = tm02*rm*tm2*[sH; vpa(subs(obj.dv)); ones(size(sV))];

                rm = [cos(obj.ba) -sin(obj.ba) 0; sin(obj.ba) cos(obj.ba) 0; 0 0 1]; % cc
                tr0 = [1 0 -sH(1); 0 1 -sV(1); 0 0 1];
                tr1 = [1 0 sH(1); 0 1 sV(1); 0 0 1];

                ta = [sH+vpa(subs(obj.dh)); sV+vpa(subs(obj.dv)); ones(size(sH))];
                tra = tr1*rm*tr0*ta;
                Hta = tra(1, 1:end);
                Vta = tra(2, 1:end);

%                 Hta = sH+Hta(2, 1:end);
%                 Vta = sV+Vta(2, 1:end);            
                
            else
                throw(MException('MATLAB:InvalidInputType','choose either <sum> or <linspace>'))
            end
        end
    end
end
