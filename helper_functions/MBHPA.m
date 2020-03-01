function [D,X,err,gerr,ratio] = ACCPlain_DL(params,varargin)


global CODE_SPARSITY CODE_ERROR codemode
global MEM_LOW MEM_NORMAL MEM_HIGH memusage
global ompfunc ompparams

CODE_SPARSITY = 1;
CODE_ERROR = 2;

MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;


%%%%% parse input parameters %%%%%


data = params.data;

ompparams = {'checkdict','off'};

% coding mode %

if (isfield(params,'codemode'))
    switch lower(params.codemode)
        case 'sparsity'
            codemode = CODE_SPARSITY;
            thresh = params.Tdata;
        case 'error'
            codemode = CODE_ERROR;
            thresh = params.Edata;
        otherwise
            error('Invalid coding mode specified');
    end
elseif (isfield(params,'Tdata'))
    codemode = CODE_SPARSITY;
    thresh = params.Tdata;
elseif (isfield(params,'Edata'))
    codemode = CODE_ERROR;
    thresh = params.Edata;
    
else
    error('Data sparse-coding target not specified');
end
ctr=0;
if (isfield(params,'trud'))
    truD=params.trud;
    ctr=1;
end
% max number of atoms %

if (codemode==CODE_ERROR && isfield(params,'maxatoms'))
    ompparams{end+1} = 'maxatoms';
    ompparams{end+1} = params.maxatoms;
end


% memory usage %

if (isfield(params,'memusage'))
    switch lower(params.memusage)
        case 'low'
            memusage = MEM_LOW;
        case 'normal'
            memusage = MEM_NORMAL;
        case 'high'
            memusage = MEM_HIGH;
        otherwise
            error('Invalid memory usage mode');
    end
else
    memusage = MEM_NORMAL;
end


% iteration count %

if (isfield(params,'iternum'))
    iternum = params.iternum;
else
    iternum = 10;
end


% omp function %

if (codemode == CODE_SPARSITY)
    ompfunc = @omp;
else
    ompfunc = @omp2;
end


% status messages %

printiter = 0;
printreplaced = 0;
printerr = 0;
printgerr = 0;

verbose = 't';
msgdelta = -1;

for i = 1:length(varargin)
    if (ischar(varargin{i}))
        verbose = varargin{i};
    elseif (isnumeric(varargin{i}))
        msgdelta = varargin{i};
    else
        error('Invalid call syntax');
    end
end

for i = 1:length(verbose)
    switch lower(verbose(i))
        case 'i'
            printiter = 1;
        case 'r'
            printiter = 1;
            printreplaced = 1;
        case 't'
            printiter = 1;
            printerr = 1;
            if (isfield(params,'testdata'))
                printgerr = 1;
            end
    end
end

if (msgdelta<=0 || isempty(verbose))
    msgdelta = -1;
end

ompparams{end+1} = 'messages';
ompparams{end+1} = msgdelta;



% compute error flag %

comperr=params.comperrdata;

% validation flag %

testgen = 0;
if (isfield(params,'testdata'))
    testdata = params.testdata;
    if (nargout>=4 || printgerr)
        testgen = 1;
    end
end


% data norms %

XtX = []; XtXg = [];
if (codemode==CODE_ERROR && memusage==MEM_HIGH)
    XtX = colnorms_squared(data);
    if (testgen)
        XtXg = colnorms_squared(testdata);
    end
end


% mutual incoherence limit %

if (isfield(params,'muthresh'))
    muthresh = params.muthresh;
else
    muthresh = 0.99;
end
if (muthresh < 0)
    error('invalid muthresh value, must be non-negative');
end


% determine dictionary size %

if (isfield(params,'initdict'))
    if (any(size(params.initdict)==1) && all(iswhole(params.initdict(:))))
        dictsize = length(params.initdict);
    else
        dictsize = size(params.initdict,2);
    end
end
if (isfield(params,'dictsize'))    % this superceedes the size determined by initdict
    dictsize = params.dictsize;
end

if (size(data,2) < dictsize)
    error('Number of training signals is smaller than number of atoms to train');
end


% initialize the dictionary %

if (isfield(params,'initdict'))
    if (any(size(params.initdict)==1) && all(iswhole(params.initdict(:))))
        D = data(:,params.initdict(1:dictsize));
    else
        if (size(params.initdict,1)~=size(data,1) || size(params.initdict,2)<dictsize)
            error('Invalid initial dictionary');
        end
        D = params.initdict(:,1:dictsize);
    end
else
    data_ids = find(colnorms_squared(data) > 1e-6);   % ensure no zero data elements are chosen
    perm = randperm(length(data_ids));
    D = data(:,data_ids(perm(1:dictsize)));
end


% normalize the dictionary %

D = normc(D);

err = zeros(1,iternum);
gerr = zeros(1,iternum);

if (codemode == CODE_SPARSITY)
    errstr = 'RMSE';
else
    errstr = 'mean atomnum';
end



%%%%%%%%%%%%%%%%%  main loop  %%%%%%%%%%%%%%%%%
Y=data;
ratio=zeros(1,iternum);
rho=params.rho;
lam=params.lam;
K=dictsize;
X=zeros(K,size(Y,2));
% initer=params.initer;
tau=params.tau;

for iter = 1:iternum
    
    G = [];
    if (memusage >= MEM_NORMAL)
        G = D'*D;
    end
    
    mu_x=rho*svds(D,1)^2;
    diff=inf;

    while diff > tau
        Xo=X;
        X=X-(1/mu_x)*D'*(D*X-Y);
        X=X.*sign(max(abs(X)-sqrt(2*lam/mu_x),0));
        diff=norm(X-Xo,'fro')/norm(Xo,'fro');

%         X=sign(X).*max(abs(X)-2*lam/mu_x,0);
    end

    %%%%%  dictionary update  %%%%%
    
    E=Y-D*X;
    suppX=abs(X)>1e-5;
    
    for j = 1:K
        d=D(:,j);
        x=X(j,:);
        E=E+d*x;
        d=E*x';
        d=d/norm(d);
        D(:,j)=d;
        E=E-d*x;
    end
    
    for k=1:K
        suppx=suppX(k,:);
        Er=E(:,suppx);
        X(k,suppx)=D(:,k)'*Er;
    end
    
    %%%%%  compute error  %%%%%
    
    
    if (testgen)
        if (memusage >= MEM_NORMAL)
            G = D'*D;
        end
        GammaG = sparsecode(testdata,D,XtXg,G,thresh);
        gerr(iter) = compute_err(D,GammaG,testdata);
    end
    
    
    %%%%%  clear dictionary  %%%%%
    
    %   D = I_clearDictionary(D,X,Y);
    
    
    %%%%%  print info  %%%%%
    
%     info = sprintf('Iteration %d / %d complete', iter, iternum);
%     if (printerr)
%         info = sprintf('%s, %s = %.4g', info, errstr, err(iter));
%     end
%     if (printgerr)
%         info = sprintf('%s, test %s = %.4g', info, errstr, gerr(iter));
%     end
%     if (printreplaced)
%         info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
%     end
%     
%     if (printiter)
%         disp(info);
%         if (msgdelta>0), disp(' '); end
%     end
%     
    if ctr
        [rat,~] = I_findDistanseBetweenDictionaries(truD,D);
        ratio(iter)=rat;
    end
    
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             sparsecode               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Gamma = sparsecode(data,D,XtX,G,thresh)

global CODE_SPARSITY codemode
global MEM_HIGH memusage
global ompfunc ompparams

if (memusage < MEM_HIGH)
    Gamma = ompfunc(D,data,G,thresh,ompparams{:});
    
else  % memusage is high
    
    if (codemode == CODE_SPARSITY)
        Gamma = ompfunc(D'*data,G,thresh,ompparams{:});
        
    else
        Gamma = ompfunc(D'*data,XtX,G,thresh,ompparams{:});
    end
    
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             compute_err              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err = compute_err(D,Gamma,data)

global CODE_SPARSITY codemode

if (codemode == CODE_SPARSITY)
    err = sqrt(sum(reperror2(data,D,Gamma))/numel(data));
else
    err = nnz(Gamma)/size(data,2);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           cleardict                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)
T2 = 0.99;
T1 = 3;
K=size(Dictionary,2);
Er=sum((Data-Dictionary*CoefMatrix).^2,1); % remove identical atoms
G=Dictionary'*Dictionary; G = G-diag(diag(G));
for jj=1:1:K,
    if max(G(jj,:))>T2 | length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 ,
        [val,pos]=max(Er);
        Er(pos(1))=0;
        Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));
        G=Dictionary'*Dictionary; G = G-diag(diag(G));
    end;
end;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            misc functions            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err2 = reperror2(X,D,Gamma)

% compute in blocks to conserve memory
err2 = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
    blockids = i : min(i+blocksize-1,size(X,2));
    err2(blockids) = sum((X(:,blockids) - D*Gamma(:,blockids)).^2);
end

end


function Y = colnorms_squared(X)

% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
    blockids = i : min(i+blocksize-1,size(X,2));
    Y(blockids) = sum(X(:,blockids).^2);
end

end
