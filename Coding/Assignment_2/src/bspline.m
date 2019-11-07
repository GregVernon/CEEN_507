classdef bspline
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        degree
        continuityVector
        numcells
        knotVector
        uniqueKnotsVector
        referenceNodes
        nodes % Control Points
        elementConnectivity
        basis
        decomposition
    end
    
    methods
        function obj = bspline(splineSpace)
            assert(length(splineSpace.uniqueKnotsVector) == length(splineSpace.continuityVector))
            
            obj = obj.splineSpace_to_knotVector(splineSpace);
            obj = obj.constructBasis();
            obj = obj.bezierExtraction("Scott");
        end
        
        function obj = splineSpace_to_knotVector(obj,splineSpace)
            obj.degree = splineSpace.degree;
            obj.uniqueKnotsVector = splineSpace.uniqueKnotsVector;
            obj.numcells = length(obj.uniqueKnotsVector)-1;
            obj.continuityVector = splineSpace.continuityVector;
            knotVector = cell(1,length(obj.uniqueKnotsVector));
            for ii = 1:length(obj.uniqueKnotsVector)
                knotVector{ii} = repmat(obj.uniqueKnotsVector(ii),1,obj.degree-obj.continuityVector(ii));
            end
            obj.knotVector = cell2mat(knotVector);
        end
        
        function obj = constructBasis(obj)
            obj.basis.variate = sym("x","real");
            
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
        
        function [obj,bez,T] = bezierExtraction(obj,method)
            if method == "Hughes"
                newContinuity = [-1 -1*ones(1,length(obj.uniqueKnotsVector)-2) -1];
            elseif method == "Scott"
                newContinuity = [-1   zeros(1,length(obj.uniqueKnotsVector)-2) -1];
            end
            newKnotVector = repelem(obj.uniqueKnotsVector,[obj.degree-newContinuity]);
            
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
            splineSpace.uniqueKnotsVector = obj.uniqueKnotsVector;
            splineSpace.continuityVector = newContinuity;
            
            bez = obj.splineSpace_to_knotVector(splineSpace);
            bez = bez.constructBasis();
            
            obj.decomposition.method = method;
            obj.decomposition.spline = bez;
            obj.decomposition.globalExtractionOperator = transpose(T{end});
            obj = collectLocalExtractionOperators(obj);
        end
        
        function [obj, C] = collectLocalExtractionOperators(obj)
            N = obj.basis.functions;
            B = obj.decomposition.spline.basis.functions;
            for e = 1:obj.numcells
                bezNodes{e} = linspace(obj.uniqueKnotsVector(e),obj.uniqueKnotsVector(e+1),obj.degree+1);
            end
            obj.decomposition.spline.nodes = unique(cell2mat(bezNodes));
            
            N_conditions = cell(length(N),1);
            for ii = 1:length(N)
                N_parts = children(N(ii));
                N_conditions{ii} = N_parts(1:end-1,2);
            end
            
            isBasisSupported = cell(obj.numcells,1);
            supportedBases = cell(obj.numcells,1);
            for e = 1:obj.numcells
                elemCenter = sum(obj.uniqueKnotsVector(e:e+1))/2;
                inDomain = false(length(N_conditions),1);
                for ii = 1:length(N_conditions)
                    condition = N_conditions{ii};
                    inDomain(ii) = any(isAlways(subs(condition,symvar(condition),elemCenter)));
                end
                isBasisSupported{e} = inDomain;
                supportedBases{e} = find(inDomain);
            end
            
            M = obj.decomposition.globalExtractionOperator;
            C = cell(obj.numcells,1);
            nodes = cell(obj.numcells,1);
            for e = 1:obj.numcells
                supportedBezierBases{e} = [((e-1)*obj.degree + 1) : (e*obj.degree+1)];
                C{e} = M(supportedBases{e},supportedBezierBases{e});
                nodes{e} = inv(C{e}')*obj.decomposition.spline.nodes(supportedBezierBases{e})';
            end
            obj.decomposition.localExtractionOperator = C;
            obj.decomposition.localExtractionSupportedSplineBases = supportedBases;
            obj.decomposition.localExtractionSupportedBezierBases = supportedBezierBases;
            obj.elementConnectivity = [obj.decomposition.localExtractionSupportedSplineBases{:}];
            obj.nodes = unique([nodes{:}])';
        end
    end
end