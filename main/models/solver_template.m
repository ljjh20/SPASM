clf

% conversion factor
mm_to_m = 1/1000;

% material property constants
constants = containers.Map;
constants('youngs_mod') = 29*10^9;
constants('shear_mod') = 10*10^9;
constants('shape_factor') = 6/5;

% Force at cantilever end
forces = containers.Map();
F = 1400;
Fa = pi; 
V = F*sin(Fa); 
H = F*cos(Fa);
forces('f') = [V; H];

% STRAIGHT BEAM SECTION
s1 = containers.Map;
s1('thickness') = 5.4*mm_to_m;
s1('length') = (125)*mm_to_m;
s1('width') = 40*mm_to_m;
s1('beam_angle') = -10*(pi/180); % rad
s1('couple') = 0; % init
s1('boundary_angle') = 0;
s1('boundary_defV') = 0;
s1('boundary_defH') = 0;
straight1 = Straight(constants, s1);

% FOR A CURVED SECTION, THE  start_angle IS THE LOADED END
c2 = containers.Map();
c2('length') = 57*mm_to_m;
c2('height') = 82*mm_to_m;
c2('thickness') = 5.4*mm_to_m;
c2('width') = 40*mm_to_m;
c2('chord') = sqrt(c2('height')^2+c2('length')^2);

% FIND SWEEP ANGLE, RADIUS
% if radius is unknown
% tangent method approximately true if length>height
% c2('sweep_angle') = 2*atan(c2('height')/c2('length'));
% c2('radius') = 0.5*c2('chord')/asin(0.5*c2('sweep_angle'));

% otherwise solve with simultaneous equations:
% length = R*(cos(a1)-cos(a2))
% height = R(sin(a2)-sin(a1))
% For desired curvature across an l,h, fix R and solve for a1, a2.
% For a desired l,h, fix a1 at starting point and solve for a2, R

% if radius is known
c2('radius') = 87.4824*mm_to_m;
c2('sweep_angle') = 2*asin(0.5*c2('chord')/c2('radius'));

c2('start_angle') = 0;
c2('end_angle') = c2('sweep_angle');
c2('boundary_angle') = 0;
c2('boundary_defV') = 0;
c2('boundary_defH') = 0;
c2('couple') = 0;
curved2 = Curved(constants, c2);

% include any number of beam sections in order they appear as entries in
% sections Map with rising integer keys
sections = containers.Map();
sections('1') = straight1;
sections('2') = curved2;

% pass hashmap into stitch function
[x_base, y_base, x_deformed, y_deformed, moment, ang, strain_e] = stitch(sections, forces);

% plotting outputs (optional)
title('deflection')
plot(x_base, y_base)
hold on
plot(x_deformed, y_deformed)
xlabel('m')
ylabel('m')
legend('static', 'deformed', location='southeastoutside')
axis('equal')

hold off
plot(x_base, moment)
title('moment') 
plot(x_base, ang)
title('angle')
plot(x_base, strain_e)
title('energy')

function [x_base, y_base, x_deformed, y_deformed, M, ang, strain_e] = stitch(sections, forces)
    f = forces('f'); % is in the form [V; H]
    
    v_offset = 0; % buffer to be overwritten for each loop
    h_offset = 0;
    a_offset = 0;
    
    M = [];
    strain_e = [];
    ang = [];
    x_base = []; % build undeformed x state of beam for plotting
    y_base = []; % build undeformed y state of beam for plotting

    x_offset = 0;
    xdef_offset = 0;
    y_offset = 0;
    ydef_offset = 0;

    x_deformed = []; % absolute x position of deformed beam
    y_deformed = []; % absolute y position of deformed beam

    angle_offset = 0;
    energy_offset= 0;

    for s = 1:size(sections)
        sec = sections(string(s)); % key in Map indicates order of sections
        height = 0; % sum as we iterate through beam sections to find moment
        length = 0;
        % sum remaining distance to point of load for moment calcualtion
        for s_other = s+1:size(sections)
            height = height + sections(string(s_other)).height;
            length = length + sections(string(s_other)).length;
        end
        if isprop(sec, "start_angle")
            if (sec.start_angle > sec.end_angle) 
                moment = cross([-height -length 0], [-f(1); -f(2); 0]); % +length
            else
                moment = cross([-height -length 0], [f(1); f(2); 0]); % +length
            end
        else
            moment = cross([-height -length 0], [f(1); f(2); 0]); % +length
        end

        sec.couple = moment(3);
        sec.ba = a_offset;
        sec.bV = v_offset;
        sec.bH = h_offset;
  
        [h_base, v_base, h_def, v_def, m, angle, energy] = sec.def(forces, 'linspace');
        h_base = h_base - h_base(1);
        v_base = v_base - v_base(1);
        h_def = h_def - h_def(1);
        v_def = v_def - v_def(1);

        M = [M m];
        strain_e = [strain_e energy+energy_offset];
        ang = [ang angle+angle_offset];

        x_base = [x_base h_base+x_offset];
        y_base = [y_base v_base+y_offset];

        x_deformed = [x_deformed h_def+xdef_offset]; % absolute x position of deformed beam
        y_deformed = [y_deformed v_def+ydef_offset]; % absolute y position of deformed beam

        x_offset = x_offset + h_base(end);
        y_offset = y_offset + v_base(end);
        xdef_offset = xdef_offset + h_def(end);
        ydef_offset = ydef_offset + v_def(end);

        angle_offset = angle(end);
        energy_offset = strain_e(end);

        a_offset = a_offset + angle(end);   
        v_offset = v_def(end) - v_base(end);
        h_offset = h_def(end) - h_base(end);

    end
end
