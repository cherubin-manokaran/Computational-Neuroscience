function [new_vector] = expandbin(old_vector, old_dt, new_dt)
    old_length = length(old_vector);
    
    scaling_ratio = new_dt / old_dt;
    
    new_length = old_length / scaling_ratio;
    
    new_vector = zeros(1, new_length);
    
    sum = 0;
    j = 1;
    for i = 1:old_length
        sum = sum + old_vector(i);
        if (mod(i, scaling_ratio) == 0)
            new_vector(j) = sum / scaling_ratio;
            j = j + 1;
            sum = 0;
        end
    end
end