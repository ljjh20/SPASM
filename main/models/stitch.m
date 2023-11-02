function [x_base, y_base, x_deformed, y_deformed, M, ang, strain_e, stiffness] = stitch(sections, forces, out_type)
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

        if (matches('curved', sec.shape))
            if (sec.start_angle > sec.end_angle) % || (obj.start_angle < 0 || obj.end_angle < 0)
                moment = cross([-height -length 0], [-f(1); -f(2); 0]); % +length
            else
                moment = cross([-height -length 0], [f(1); f(2); 0]); % +length
            end
        elseif (matches('straight', sec.shape))
            if (sec.length < 0) % || (obj.start_angle < 0 || obj.end_angle < 0)
                moment = cross([-height -length 0], [-f(1); -f(2); 0]); % +length
            else
                moment = cross([-height -length 0], [f(1); f(2); 0]); % +length
            end
        end

        sec.couple = moment(3);
        sec.ba = a_offset;
%         sec.ba = 0;
%         sec.bV = v_offset;
%         sec.bH = h_offset;

        if(matches('linspace', out_type))
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
        elseif(matches('sum', out_type))
            [h_base, v_base, h_def, v_def, m, angle, energy] = sec.def(forces, 'sum');
    
            M = [M m];
            strain_e = [strain_e energy];
            ang = [ang angle];
    
            x_base = [x_base h_base];
            y_base = [y_base v_base];
    
            x_deformed = [x_deformed h_def]; % absolute x position of deformed beam
            y_deformed = [y_deformed v_def]; % absolute y position of deformed beam
   
            a_offset = a_offset + angle;  
%             v_offset = v_def(end) - v_base(end);
%             h_offset = h_def(end) - h_base(end);
        else
            throw(MException('MATLAB:InvalidInputType','choose either <sum> or <linspace>'))
        end
    end 

    if(matches('sum', out_type))
        x_base = sum(x_base);
        y_base = sum(y_base);
        x_deformed = sum(x_deformed);
        y_deformed = sum(y_deformed);
        M = M(1); % couple is passed through internally
        ang = sum(ang);
        strain_e = sum(strain_e);
    end
    stiffness = M(1)/ang(end);
end