classdef linear_eq
    properties
        vector         % element of Xi_vector or small_Xi_vector
        constant       % array with the number of nodes
    end
    methods
        function w = linear_eqs(coefs, const)
           % function w = linear_eqs(coefs, const)
           if ~isa(coefs,'Xi_vector') && ~isa(coefs,'small_Xi_vector') 
               error('The coefficients have to be of the prescribed form')
           end
           if ~isa(const, 'float')
               error('The constant has to be a single floating point element')
           end
           w.vector = coefs;
           w.constant = const;
        end
        
        function c = eval(lin_eq, vec)
            % function c = eval(lin_eq, vec)
            % evaluate < vec, lin_eq > + c

            if ~isa(lin_eq,'linear_eq')
                temp = lin_eq;
                lin_eq = vec;
                vec = temp;
            end
            if ~compatible(lin_eq.vector, vec)
                error('The coefficients of the linear equation and the vector are incmpatible')
            end
            
            c = lin_prod(vec, lin_eq.vector) + lin_eq.const;
        end
        
        function d = der(lin_eq)
            d = Xi2vec(lin_eq.vector);
        end
    end
end
