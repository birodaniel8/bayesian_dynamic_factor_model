classdef BDFM < handle
    %BDFM - Bayesian Dynamic Factor Model class
    %   Detailed explanation goes here

    properties
        y, m, T, k
    end

    methods
        function obj = BDFM(y, k)
            %BDFM Construct an instance of BDFM
            %   Detailed explanation goes here
            obj.y = y;
            [obj.T, obj.m] = size(y);
            obj.k = k;
        end
    end
end
