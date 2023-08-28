function varargout = liouvilleGuiM(varargin)
% LIOUVILLEGUIM MATLAB code for liouvilleGuiM.fig
%      LIOUVILLEGUIM, by itself, creates a new LIOUVILLEGUIM or raises the existing
%      singleton*.
%
%      H = LIOUVILLEGUIM returns the handle to a new LIOUVILLEGUIM or the handle to
%      the existing singleton*.
%
%      LIOUVILLEGUIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LIOUVILLEGUIM.M with the given input arguments.
%
%      LIOUVILLEGUIM('Property','Value',...) creates a new LIOUVILLEGUIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before liouvilleGuiM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to liouvilleGuiM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help liouvilleGuiM

% Last Modified by GUIDE v2.5 21-Oct-2021 14:44:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @liouvilleGuiM_OpeningFcn, ...
                   'gui_OutputFcn',  @liouvilleGuiM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before liouvilleGuiM is made visible.
function liouvilleGuiM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to liouvilleGuiM (see VARARGIN)


% Choose default command line output for liouvilleGuiM
handles.output = hObject;
draw_axes(handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes liouvilleGuiM wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end

% --- Outputs from this function are returned to the command line.
function varargout = liouvilleGuiM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

function RhoDimension_Callback(hObject, eventdata, handles)
% hObject    handle to RhoDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: get(hObject,'String') returns contents of RhoDimension as text
%        str2double(get(hObject,'String')) returns contents of RhoDimension as a double
end

% --- Executes during object creation, after setting all properties.
function RhoDimension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RhoDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function NumQubit_Callback(hObject, eventdata, handles)
% hObject    handle to NumQubit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumQubit as text
%        str2double(get(hObject,'String')) returns contents of NumQubit as a double


% --- Executes during object creation, after setting all properties.
end
function NumQubit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumQubit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function listOf_H_Callback(hObject, eventdata, handles)
% hObject    handle to listOf_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of listOf_H as text
%        str2double(get(hObject,'String')) returns contents of listOf_H as a double
end

% --- Executes during object creation, after setting all properties.
function listOf_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listOf_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function MaxErrorPower_Callback(hObject, eventdata, handles)
% hObject    handle to MaxErrorPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxErrorPower as text
%        str2double(get(hObject,'String')) returns contents of MaxErrorPower as a double

end
% --- Executes during object creation, after setting all properties.
function MaxErrorPower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxErrorPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function InitRho_Callback(hObject, eventdata, handles)
% hObject    handle to InitRho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of InitRho as text
%        str2double(get(hObject,'String')) returns contents of InitRho as a double
end

% --- Executes during object creation, after setting all properties.
function InitRho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InitRho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

end
% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function DecoherenceMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to DecoherenceMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DecoherenceMatrix as text
%        str2double(get(hObject,'String')) returns contents of DecoherenceMatrix as a double

end
% --- Executes during object creation, after setting all properties.
function DecoherenceMatrix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DecoherenceMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pl1.
function pl1_Callback(hObject, eventdata, handles)
% hObject    handle to pl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recall(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of pl1
end


% --- Executes on button press in pl2.
function pl2_Callback(hObject, eventdata, handles)
% hObject    handle to pl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recall(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of pl2
end

% --- Executes on button press in pl3.
function pl3_Callback(hObject, eventdata, handles)
% hObject    handle to pl3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recall(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of pl3
end
% --- Executes on button press in zero_H.
function zero_H_Callback(hObject, eventdata, handles)
% hObject    handle to zero_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recall(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of zero_H
end

% --- Executes on button press in Ok.
function Ok_Callback(hObject, eventdata, handles)
% hObject    handle to Ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recall(hObject, eventdata, handles)
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function rho_total = initialState(initMatrix, numQubits)
    rho_total = initMatrix;
    for i=2: numQubits
        rho_total = kron(rho_total,initMatrix);
    end
end
% -----------------------------------------------------
% Liouvillian .  --------------------------------------
function L = Liouvillian(numQubits, decoheranceM)

    A_t = cell(1,numQubits);
    for i = 1: numQubits
        A_t{i} = eye(2);
    end
    L = zeros(4^numQubits,4^numQubits);
    for i=1:numQubits
        A_t{i} = decoheranceM;
        A = tensorPro(A_t, numQubits);
        A_dag = ctranspose(A);
        I = eye(2^numQubits);
        L = L + kron(A, A_dag) - 0.5 * (kron(I, A_dag * A ) + kron(A_dag * A, I));
        A_t{i} = eye(2);
    end
end

function T = tensorPro(A, numQubits)
    T = A{1};
    for i=2:numQubits
        T = kron(T,A{i});
    end
end
% -----------------------------------------------------
% Display Graph  --------------------------------------
function s = startIndex(cir)
    global db
    s = db.max_error_power*(cir-1) + 1;

end
% -----------------------------------------------------
% Display Graph  --------------------------------------
function disp_graph(x, yA)
    plot(x,yA')
end
% -----------------------------------------------------
% Generating random hamiltonian  ----------------------
function Hl = R_Hamiltonian(dim)   
    global db
    Ones = ones(dim);
    C_rel = 2 * pi * rand(dim,dim)- pi * Ones;
    C_img = 2 * pi * rand(dim,dim)- pi * Ones;
    C = C_rel + C_img * j;
    C_t = ctranspose(C);
    H = C + C_t; 
    if db.zeroH
        H = zeros(dim,dim);
    end
    I = eye(dim);
    Hl = kron(H,I) - kron(I,transpose(H));

end
% -----------------------------------------------------
% kraus operator---------------------------------------
function [k,kp] = kraus(error_power,liovilian, ham_dim,m,n)     
    global db
    R_H = R_Hamiltonian(ham_dim);
    k = expm(j * R_H + (liovilian * error_power));
    kp = expm((-j * R_H + (liovilian * error_power)));
    db.eignValueVec(((n-1)*(db.max_error_power))+m) = (min(eig(R_H)) + max(eig(R_H)))/2;
%     db.eignValueVecMax(((n-1)*(db.max_error_power))+m) = min(eig(R_H));

end
% -----------------------------------------------------
% Rho_final operate on rho0 as victor   --  A  --------
function rho_f = rho_f2(num_circuits, rho0, dim, error_power, liovilian,m,n) 

    [k,kp] = kraus(error_power,liovilian, dim,m,n);
    K = k * kp;
    for i=2  : num_circuits
        [k,kp] = kraus(error_power,liovilian, dim,m,n);
        K = k * K * kp;
    end
    rho_f = K * rho0';
end
% -----------------------------------------------------
% Rho_final operate on rho0 as victor   --  B  --------
function rho_f = rho_f1(num_circuits, rho0, dim, error_power, liovilian,m,n) 
    
    [,k_] = kraus(error_power,liovilian, dim,m,n);
    for i= 2 : num_circuits
        [,k] = kraus(error_power,liovilian, dim,m,n);
        k_ = k * k_;
    end
    rho_f = k_ * rho0';

end
% -----------------------------------------------------
% trace(rho-f^2)  -------------------------------------
function pur = purity(rho_f , dimRho) 

    rho_f_matrix = reshape(rho_f, dimRho, dimRho);
    pur = trace(rho_f_matrix * rho_f_matrix);
    
end
function list = getList(str0)    
    str1 = extractBefore(str0,strlength(str0));
    str = extractAfter(str1,1);
    
    my =  split(str," ");
    list = zeros(1,length(my));
    for i = 1 : length(list)
        list(i) = str2num(my{i});
    end
end
% -----------------------------------------------------
% convert string matrix to matrix  --------------------
function M = matrix(str0, dim1,dim2)
    str1 = extractBefore(str0,strlength(str0));
    str = extractAfter(str1,1);
    my =  split(str,";");
    M = zeros(dim1,dim2);
    for i = 1 : length(my)
        ele = split(my(i)," ");
        for j = 1 : length(ele)
            M(j+ (i-1)*dim2) =  str2double(regexp(ele{j},'[+-]?\d+','match'));
        end
    end
    if dim1 == dim2
        M = transpose(M);
    end
end

function recall(hObject, eventdata, handles)
    global db
    
    db.randomHamilList =  getList(get(handles.listOf_H,'String'));
    db.num_circuits =  length(db.randomHamilList);
    db.num_qubit = str2num(get(handles.NumQubit,'String'));
    db.rho0 = get(handles.InitRho,'String');
    db.A_ = get(handles.DecoherenceMatrix,'String');
    db.max_error_power = str2num(get(handles.MaxErrorPower,'String'));
    db.zeroH = get(handles.zero_H, 'Value');
    
    h_dim = 2^db.num_qubit;
    dim_rho = 2^db.num_qubit;
    
    rho_0_2D = (0.5) * matrix(db.rho0,2,2);
    rho_0 = initialState(rho_0_2D, db.num_qubit);% as matrix (tensor product of the initial state rho1 tensor rho1)
    rho_0 = reshape(rho_0',4^db.num_qubit,1);% as horizental vector
    rho_0 = rho_0';% as vertical vector

    A =  0.0001 * matrix(db.A_,2,2);
    L = Liouvillian(db.num_qubit, A);

    liovilian =   L;
    db.power_error_vec = zeros(1,db.max_error_power);
    db.ratio_vecA = zeros(db.max_error_power, db.num_circuits);% A - 2B
    db.ratio_vecB = zeros(db.max_error_power, db.num_circuits);%(A/2B)-1
    db.ratio_vecC = zeros(db.max_error_power, db.num_circuits);%Puirty <KKp*rho|KKp*rho>
    db.ratio_vecD = zeros(db.max_error_power, db.num_circuits);%Puirty <K*rho|K*rho>
    db.power_error_vec = 1:db.max_error_power;
    db.eignValueVec = zeros(db.max_error_power, db.num_circuits);% eignvalue
    db.traceOfRhoF = zeros(db.max_error_power, db.num_circuits);% eignvalue

    for j=1 : db.num_circuits

        for i=1 : db.max_error_power
            rhof1 = rho_f1(db.randomHamilList(j), rho_0, h_dim, i, liovilian,i,j) ;
            rhof2 = rho_f2(db.randomHamilList(j), rho_0, h_dim, i, liovilian,i,j) ;
            traceOfRho_F = trace(reshape(rhof1, dim_rho, dim_rho));
            pur1 = purity(rhof1, dim_rho);
            pur2 = purity(rhof2, dim_rho);
            
            ratio1 = (1-pur2) -  2 * (1-pur1);
            ratio2 = 1- ((1-pur2) / 2 * (1-pur1));

            db.ratio_vecA(((j-1)*(db.max_error_power))+i) = ratio1;
            db.ratio_vecB(((j-1)*(db.max_error_power))+i) = ratio2;
            db.ratio_vecC(((j-1)*(db.max_error_power))+i) = pur2;
            db.ratio_vecD(((j-1)*(db.max_error_power))+i) = pur1;
            db.traceOfRhoF(((j-1)*(db.max_error_power))+i) = traceOfRho_F;
        end
    end
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------

    s1 = startIndex(1);
    s2 = startIndex(2);
    s3 = startIndex(3);
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
    cla(handles.axes1)
    axes(handles.axes1)
    if get(handles.pl1,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecA(s1:s1+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl2,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecA(s2:s2+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl3,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecA(s3:s3+db.max_error_power-1))
    end

% -------------- Graph  2   -----------------------------------------------

    cla(handles.axes2)
    axes(handles.axes2)
    if get(handles.pl1,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecB(s1:s1+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl2,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecB(s2:s2+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl3,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecB(s3:s3+db.max_error_power-1))
    end

% -------------- Graph  3   -----------------------------------------------
  
    cla(handles.axes3)
    axes(handles.axes3)
    if get(handles.pl1,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecC(s1:s1+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl2,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecC(s2:s2+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl3,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecC(s3:s3+db.max_error_power-1))
    end

% -------------- Graph  4   -----------------------------------------------

    cla(handles.axes4)
    axes(handles.axes4)
    if get(handles.pl1,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecD(s1:s1+db.max_error_power-1))
    end 
    hold on
    if get(handles.pl2,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecD(s2:s2+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl3,'Value')
        disp_graph(db.power_error_vec,db.ratio_vecD(s3:s3+db.max_error_power-1))
    end

% -------------- Graph  5   -----------------------------------------------

    cla(handles.axes5)
    axes(handles.axes5)
    if get(handles.pl1,'Value')
        disp_graph(db.power_error_vec,db.eignValueVec(s1:s1+db.max_error_power-1))
    end 
    hold on
    if get(handles.pl2,'Value')
        disp_graph(db.power_error_vec,db.eignValueVec(s2:s2+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl3,'Value')
        disp_graph(db.power_error_vec,db.eignValueVec(s3:s3+db.max_error_power-1))
    end
% -------------- Graph  6   -----------------------------------------------

    cla(handles.axes6)
    axes(handles.axes6)
    if get(handles.pl1,'Value')
        disp_graph(db.power_error_vec,db.traceOfRhoF(s1:s1+db.max_error_power-1))
    end 
    hold on
    if get(handles.pl2,'Value')
        disp_graph(db.power_error_vec,db.traceOfRhoF(s2:s2+db.max_error_power-1))
    end
    hold on 
    if get(handles.pl3,'Value')
        disp_graph(db.power_error_vec,db.traceOfRhoF(s3:s3+db.max_error_power-1))
    end
    
    draw_axes(handles)

end


function draw_axes(handles)

    axes(handles.axes1)
    title('        (1-<kkp\rho0|kkp\rho0>) - 2(1-<k\rho0|k\rho0>))  VS Error Power','LineWidth',14)
    xlabel('ERROR POWER')
    ylabel('(1-<kkp\rho0|kkp\rho0>) - 2(1-<k\rho0|k\rho0>)) ')
    if get(handles.zero_H, 'Value')
        legend('zero Hmailtonian','zero Hmailtonian','zero Hmailtonian')
    else
        legend('One Random circuit','Two Random circuits','Three Random circuits')
    end

    axes(handles.axes2)
    title(' 1-[(1-<kkp\rho0|kkp\rho0> / 2(1-<k\rho0|k\rho0>))] VS Error Power')
    xlabel('ERROR POWER')
    ylabel('1-(1-<kkp\rho0|kkp\rho0> / 2(1-<krho0|k\rho0>))')
    if get(handles.zero_H, 'Value')
        legend('zero Hmailtonian','zero Hmailtonian','zero Hmailtonian')
    else
        legend('One Random circuit','Two Random circuits','Three Random circuits')
    end
    
    axes(handles.axes3)
    title('<kkp\rho0|kkp\rho0> Vs Error Power')
    xlabel('ERROR POWER')
    ylabel('<kkp\rho0|kkp\rho0>')
    if get(handles.zero_H, 'Value')
        legend('zero Hmailtonian','zero Hmailtonian','zero Hmailtonian')
    else
        legend('One Random circuit','Two Random circuits','Three Random circuits')
    end
    
    axes(handles.axes4)
    title('<k\rho0|k\rho0> Vs Error Power')
    xlabel('ERROR POWER')
    ylabel('<k\rho0|k\rho0>')
    if get(handles.zero_H, 'Value')
        legend('zero Hmailtonian','zero Hmailtonian','zero Hmailtonian')
    else
        legend('One Random circuit','Two Random circuits','Three Random circuits')
    end
    
    axes(handles.axes5)
    title('Hamiltonian eignvalue mean Vs Error Power')
%     xlabel('ERROR POWER')
    ylabel('Hamiltonian eignvalue mean')
    if get(handles.zero_H, 'Value')
        legend('zero Hmailtonian','zero Hmailtonian','zero Hmailtonian')
    else
        legend('One Random circuit','Two Random circuits','Three Random circuits')
    end

    axes(handles.axes6)
    title('Trace preserving')
    xlabel('ERROR POWER')
    ylabel('Trace of final state')
    if get(handles.zero_H, 'Value')
        legend('zero Hmailtonian','zero Hmailtonian','zero Hmailtonian')
    else
        legend('One Random circuit','Two Random circuits','Three Random circuits')
    end
end
