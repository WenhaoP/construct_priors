function [mask] = makeSeaMask(loc, tv, border)
    % Check if mesh node is inside Norway
    mask1 = inpolygon(loc(:,1), loc(:,2), border(:,1), border(:,2)); 
    
    % Expand one level
    mask = mask1;
    for i = 1:size(tv, 1)
        if(mask1(tv(i,1)) == 1)
            mask(tv(i,2)) = 1;
            mask(tv(i,3)) = 1;
        end
        if(mask1(tv(i,2)) == 1)
            mask(tv(i,1)) = 1;
            mask(tv(i,3)) = 1;
        end
        if(mask1(tv(i,3)) == 1)
            mask(tv(i,1)) = 1;
            mask(tv(i,2)) = 1;
        end        
    end
    
    % Expand one level
    mask1 = mask;
    mask = mask1;
    for i = 1:size(tv, 1)
        if(mask1(tv(i,1)) == 1)
            mask(tv(i,2)) = 1;
            mask(tv(i,3)) = 1;
        end
        if(mask1(tv(i,2)) == 1)
            mask(tv(i,1)) = 1;
            mask(tv(i,3)) = 1;
        end
        if(mask1(tv(i,3)) == 1)
            mask(tv(i,1)) = 1;
            mask(tv(i,2)) = 1;
        end        
    end
end
