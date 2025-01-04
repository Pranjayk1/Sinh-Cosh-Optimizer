%% Copyright (c) 2020, Mehdi Ghasri
%% To cooperate in articles, send an email to the following address
% (with Subject=CO Article):
% Email: Eng.mehdighasri@gmail.com


%% Copyright (c) 2019, Mehdi Ghasri
% All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of IUPUI nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [lb,ub,dim,fobj] = CEC2017(F)


switch F
    %CEC2017 - 1 
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 2
    case 'F2'
        fobj = @F2;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 3
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 4
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 5
    case 'F5'
        fobj = @F5;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 6
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 7
    case 'F7'
        fobj = @F7;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 8
    case 'F8'
        fobj = @F8;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 9
    case 'F9'
        fobj = @F9;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 10
    case 'F10'
        fobj = @F10;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 11
    case 'F11'
        fobj = @F11;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 12
    case 'F12'
        fobj = @F12;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 13
    case 'F13'
        fobj = @F13;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 14
    case 'F14'
        fobj = @F14;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 15
    case 'F15'
        fobj = @F15;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 16
    case 'F16'
        fobj = @F16;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 17
    case 'F17'
        fobj = @F17;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 18
    case 'F18'
        fobj = @F18;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 19
    case 'F19'
        fobj = @F19;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 20
    case 'F20'
        fobj = @F20;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 21
    case 'F21'
        fobj = @F21;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 22
    case 'F22'
        fobj = @F22;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 23
    case 'F23'
        fobj = @F23;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 24
    case 'F24'
        fobj = @F24;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 25
    case 'F25'
        fobj = @F25;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 26
    case 'F26'
        fobj = @F26;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 27
    case 'F27'
        fobj = @F27;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 28
    case 'F28'
        fobj = @F28;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 29
    case 'F29'
        fobj = @F29;
        lb=-100;
        ub=100;
        dim=10;
        
    %CEC2017 - 30
    case 'F30'
        fobj = @F30;
        lb=-100;
        ub=100;
        dim=10;
end

end

% F1
function o = F1(x) 
    o = cec17_func(x',1);





    
end

% F2
function o = F2(x) 
    o = cec17_func(x',2);
end

% F3
function o = F3(x) 
    o = cec17_func(x',3);
end

% F4
function o = F4(x) 
    o = cec17_func(x',4);
end

% F5
function o = F5(x) 
    o = cec17_func(x',5);
end

% F6
function o = F6(x) 
    o = cec17_func(x',6);
end

% F7
function o = F7(x) 
    o = cec17_func(x',7);
end

% F8
function o = F8(x) 
    o = cec17_func(x',8);
end

% F9
function o = F9(x) 
    o = cec17_func(x',9);
end

% F10
function o = F10(x) 
    o = cec17_func(x',10);
end


% F11
function o = F11(x) 
    o = cec17_func(x',11);
end

% F12
function o = F12(x) 
    o = cec17_func(x',12);
end


% F13
function o = F13(x) 
    o = cec17_func(x',13);
end

% F14
function o = F14(x) 
    o = cec17_func(x',14);
end

% F15
function o = F15(x) 
    o = cec17_func(x',15);
end

% F16
function o = F16(x) 
    o = cec17_func(x',16);
end


% F17
function o = F17(x) 
    o = cec17_func(x',17);
end

% F18
function o = F18(x) 
    o = cec17_func(x',18);
end

% F19
function o = F19(x) 
    o = cec17_func(x',19);
end

% F20
function o = F20(x) 
    o = cec17_func(x',20);
end

% F21
function o = F21(x) 
    o = cec17_func(x',21);
end

% F22
function o = F22(x) 
    o = cec17_func(x',22);
end

% F23
function o = F23(x) 
    o = cec17_func(x',23);
end

% F24
function o = F24(x) 
    o = cec17_func(x',24);
end


% F25
function o = F25(x) 
    o = cec17_func(x',25);
end

% F26
function o = F26(x) 
    o = cec17_func(x',26);
end

% F27
function o = F27(x) 
    o = cec17_func(x',27);
end

% F28
function o = F28(x) 
    o = cec17_func(x',28);
end

% F29
function o = F29(x) 
    o = cec17_func(x',29);
end

% F30
function o = F30(x) 
    o = cec17_func(x',30);
end
