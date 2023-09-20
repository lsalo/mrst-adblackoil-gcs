function f = addScanKr(f, regi, minSat)
%
% January 2020: Support for relative permeability hysteresis in the
%               nonwetting phase using Killough's (1976) model. This uses
%               the Land (1968) model to determine the trapped gas
%               saturation (Sgt) and then Killough's model to determine the
%               gas relative permeability along a scanning curve. Scanning
%               curves are reversible in this model.
%
%               It requires that the bounding imbibition curve already be
%               in the fluid structure in krX, after drainage tables. 
%               Currently, this can be obtained from assignSGFN (keyword 
%               SGFN in ECLIPSE input file). Assumes that gas is the 
%               nonwetting phase.
%        

% N imbibition regions
regi  = unique(regi);
nregi = numel(regi);
reg = regi-nregi;

% Maximum and minimum gas saturation
if isfield(f, 'sWcon')
    sgmax = 1 - f.sWcon;
    if numel(sgmax) ~= nregi  % temporary
        sgmax = ones(1, nregi)*sgmax;
    end
    sgmin = zeros(1, nregi);
    warning('Min gas saturation (primary drainage) set to 0')
elseif isfield(f, 'krPts')
    f.sgtmax   = f.krPts.g(regi,2)';
    sgmax      = f.krPts.g(regi, 3)';
    sgmin      = f.krPts.g(reg, 1)';
else
    sgmax = ones(1, nregi);
    sgmin = zeros(1, nregi);
    warning('Max gas saturation (primary drainage) set to 1.')
    warning('Min gas saturation (primary drainage) set to 0')
end
assert(numel(sgmax) == nregi);

% Add to fluid object
kri = cell(1, nregi);
for n = 1:nregi
    kriVars = {sgmax(n), sgmin(n), f.krG{n}, f.krG{n+nregi}, f.sgtmax(n)};
    kri{n}  = @(sg, sgi) scanRelPerm(sg, sgi, kriVars, minSat, f.ehystr);
end
f.krGi   = kri;

% ----------------- Compute Scanning kr curve --------------------
    function [kri, tol, minSat, sgt] = scanRelPerm(sg, sgi, kriVars, minSat, opts)
        % Compute scanning relative permeability curves.
        % sg  = current gas saturation
        % sgi = maximum gas saturation achieved in the current run.
        % __________________________________________________________
        nci = numelValue(sg);
        
        % Set tolerance for flow reversal and initialize
        assert(nci == numelValue(sgi))   % same n for both sg and sgi
        tol  = 1e-3;
        kri  = nan(nci, 1);
        AD = getSampleAD(sg, sgi);
        if isa(AD, 'GenericAD')
            kri = double2GenericAD(kri, AD);
        elseif isa(AD, 'ADI')
            kri = double2ADI(kri, AD);
        end
        
        % Inputs
        hystModel = opts{2};
        sgmx      = kriVars{1};
        sgmn      = kriVars{2};
        krG       = kriVars{3};
        krGib     = kriVars{4};
        sgtmx     = kriVars{5};
        
        % Index for cells where scanning curves are computed. We set a
        % minimum of 5% gas saturation to activate hysteresis.
        imb  = all([value(sg) + tol < value(sgi), ...
                    value(sgi) > minSat+sgmn], 2);
        
        if any(imb)
            kri(~imb) = krG(sg(~imb));
            sgiv  = sgi(imb);
            sgv   = sg(imb);
            
            if hystModel == 2           % Killough (1976) nonwetting phase.
                a     = opts{4};        % modif. param. for improved converg.
                A     = 1 + a*(sgmx - sgiv);
                C     = 1/(sgtmx - sgmn) - 1/(sgmx-sgmn);
                sgt   = sgmn + (sgiv-sgmn)./(A + C*(sgiv-sgmn));
                snorm = sgtmx + (sgv-sgt)*(sgmx-sgtmx)./(sgiv-sgt);
                v = krGib(snorm).*krG(sgiv)./krG(sgmx);
                kri(imb) = v;
                
                isAbove = kri > krG(value(sg)) & imb;
                if any(isAbove)
                    kri(isAbove) = krG(sg(isAbove));
                    % kri(isAbove) = krG(value(sg(isAbove))) - 1e-5;
                end
            else
                error(['Indicated hysteresis model in item 2 of EHYSTR' ...
                    ' keyword in the .DATA input file is not supported.'])
            end
            
            % Checks
            isBelow = all([sgv < sgt, value(kri(imb)) ~= 0], 2);
            assert(~any(isBelow));
            kriv = value(kri);
            assert(all(~isnan(kriv)));                  % no nan values
            assert(all(kriv >= 0));                     % 0 or positive values
            assert(all(kriv <= krG(value(sgmx))));      % max kr not exceeded
        else
            kri = krG(sg);                      % all cells on drainage curve
        end
    end

end