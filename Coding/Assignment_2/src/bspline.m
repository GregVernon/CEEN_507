classdef bspline
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        degree
        continuityVector
        numcells
        nodes
        knotVector
        basis
    end
    
    methods
        function obj = bspline(splineSpace)
            assert(length(splineSpace.nodes) == length(splineSpace.continuityVector))
            
            obj = obj.splineSpace_to_knotVector(splineSpace);
            obj = obj.constructBasis();
        end
        
        function obj = splineSpace_to_knotVector(obj,splineSpace)
            obj.degree = splineSpace.degree;
            obj.nodes = splineSpace.nodes;
            obj.numcells = length(obj.nodes)-1;
            obj.continuityVector = splineSpace.continuityVector;
            knotVector = cell(1,length(obj.nodes));
            for ii = 1:length(obj.nodes)
                knotVector{ii} = repmat(obj.nodes(ii),1,obj.degree-obj.continuityVector(ii));
            end
            obj.knotVector = cell2mat(knotVector);
        end
        
        function obj = constructBasis(obj)
            obj.basis.variate = sym("xi","real");
            
            x = obj.basis.variate;
            kV = obj.knotVector;
            N = cell(obj.degree+1,1);
            for p = 0:obj.degree
                N{p+1} = sym(zeros(length(obj.knotVector)-(p+1),1));
                for ii = 1:length(obj.knotVector)-(p+1)
                    if p == 0
                        if ii < length(obj.knotVector)-(p+1)
                            N{p+1}(ii) = piecewise(kV(ii) <= x < kV(ii+1),sym(1),sym(0));
                        else
                            N{p+1}(ii) = piecewise(kV(ii) <= x <= kV(ii+1),sym(1),sym(0));
                        end
                    else
                        divisor(1) = (kV(ii+p)-kV(ii));
                        if divisor(1) == 0
                            term(1) = sym(0);
                        else
                            term(1) = ((x - kV(ii))/(kV(ii+p)-kV(ii)));
                        end
                        
                        divisor(2) = (kV(ii+p+1)-kV(ii+1));
                        if divisor(2) == 0
                            term(2) = sym(0);
                        else
                            term(2) = ((kV(ii+p+1) - x)/(kV(ii+p+1)-kV(ii+1)));
                        end
                        N{p+1}(ii) = term(1)*N{p}(ii) + term(2)*N{p}(ii+1);
                    end
                end
            end
            obj.basis.functions = simplify(N{end});
        end
        
        function [obj,T] = bezierExtraction(obj,method)
            if method == "Hughes"
                newContinuity = [-1 -1*ones(1,length(obj.nodes)-2) -1];
            elseif method == "Scott"
                newContinuity = [-1 -1*ones(1,length(obj.nodes)-2) -1];
            end
            newKnotVector = repelem(obj.nodes,[obj.degree-newContinuity]);
            
            KV = obj.knotVector;
            nKV = newKnotVector;
            T = cell(obj.degree+1,1);
            for q = 0:obj.degree
                nRows = ((length(nKV)-1)-q);
                nCols = ((length(KV)-1)-q);
                T{q+1} = sym(zeros(nRows,nCols));
                for ii = 1:((length(nKV)-1)-q)
                    for jj = 1:((length(KV)-1)-q)
                        if q == 0
                            if KV(jj) <= nKV(ii) && nKV(ii) < KV(jj+1)
                                T{q+1}(ii,jj) = sym(1);
                            else
                                T{q+1}(ii,jj) = sym(0);
                            end
                        else
                            D1 = KV(jj+q) - KV(jj);
                            if D1 == 0
                                A = sym(0);
                            else
                                A = (nKV(ii+q) - KV(jj)) / D1;
                            end
                            
                            D2 = KV(jj+q+1) - KV(jj+1);
                            if D2 == 0
                                B = sym(0);
                            else
                                B = (KV(jj+q+1) - nKV(ii+q)) / D2;
                            end
                            
                            T{q+1}(ii,jj) = A * T{q}(ii,jj) + B * T{q}(ii,jj+1);
                        end
                    end
                end
            end
            splineSpace.degree = obj.degree;
            splineSpace.nodes = obj.nodes;
            splineSpace.continuityVector = newContinuity;
            
            obj = bspline(splineSpace);
        end
    end
end

