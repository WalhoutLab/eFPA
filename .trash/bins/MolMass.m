% Copyright (c) 2009, Eckart Giebler
% function MolMass.m
% The function MolMass calculates the molar weight of a substance.
% 
% function call: MM = MolMass(substance)
% 
% substance is a string of the chemical formula of s substance.
% example:	MM = MolMass('Fe2(SO4)3');
% 
% substance can also be a vector of substances opened by '[' and divided by space, comma or semicolon.
% examples:
% 			MM = MolMass('[Fe2(SO4)3 CuSO4 NaOH]');
% 			MM = MolMass('[H2SO4;H2O;P;Cl2]');
% 			MM = MolMass('[C3H5(OH)3,C3H7OH,C12H22O11,NaCl]');
% 
% To distinguish charched substances the symbols '+' and '-' can be used.
% exampels:
% 			MM = MolMass('Fe2+')  --->  MM = 55.8470		%it means one mol of Fe2+
% 			MM = MolMass('Fe3+')  --->  MM = 55.8470		%it means one mol of Fe3+
% 		but	MM = MolMass('Fe2')   --->  MM = 111.6940		%it means two moles of Fe
% 
% 
% ultimate Date:     18.12.2009
% version:           1.4
% copyright:         E. Giebler, TU Dresden, IfA (1999-2009)
% Modified By: Xuhang Li: 2020 
%               - modified for existance of decimals in the formula


function MM = MolMass(substance)

for char_index = 1:length(substance),														% first syntax check
   if isempty(findstr(substance(char_index),'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[];, +-".')),
      error_message = ['Wrong symbol ' substance(char_index) ' in the formula!']
      error(error_message)
   end%if
end%for

substance_number = 1;
single_substance = '';
char_index=0;

% ---------------------replace redundant spaces--------------------------
if   ~isempty(findstr('  ',substance))...
   | ~isempty(findstr(' ]',substance))...
   | ~isempty(findstr('[ ',substance))
  while ~isempty(findstr('  ',substance))
    substance = strrep(substance,'  ',' ');
  end % while  ~isempty(findstr('  ',substance))
    substance = strrep(substance,' ]',']');
    substance = strrep(substance,'[ ','[');
end %if findstr(...

if substance(1) == '[',																		% vector of substances
char_index = char_index + 1;
   while char_index < length(substance)
      char_index = char_index + 1;
      if ~isempty(findstr(substance(char_index),'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()+-".'))
         single_substance(end+1) = substance(char_index);
      elseif ~isempty(findstr(substance(char_index),' ,;]'))
         MM(substance_number) = MolMass(single_substance);
         single_substance = '';
         substance_number = substance_number+1;
         while char_index < length(substance) &  ~isempty(findstr(substance(char_index+1),' ,;]'))
            char_index = char_index + 1;
         end%while
      	if substance(char_index) == ']'
         	char_index = length(substance);
      	end%if
      end%if
   end%while
 %  MM(substance_number) = MolMass(single_substance);
else														                                % single substance only
   
% molar masses of the elements:
H = 1.00794; He = 4.002602; Li = 6.941; Be = 9.012182; B = 10.811; C = 12.011; N = 14.00674;
O = 15.9994; F = 18.9984032; Ne = 20.1797; Na = 22.989768; Mg = 24.305; Al = 26.981539; Si = 28.0855;
P = 30.973762; S = 32.066; Cl = 35.4527; Ar = 39.948; K = 39.0983; Ca = 40.078; Sc = 44.95591;
Ti = 47.88; V = 50.9415; Cr = 51.9961; Mn = 54.93805; Fe = 55.847; Co = 58.9332; Ni = 58.69;
Cu = 63.546; Zn = 65.39; Ga = 69.723; Ge = 72.61; As = 74.92159; Se = 78.96; Br = 79.904; Kr = 83.8;
Rb = 85.4678; Sr = 87.62; Y = 88.90585; Zr = 91.224; Nb = 92.90638; Mo = 95.94; Tc = 98.9063;
Ru = 101.07; Rh = 102.9055; Pd = 106.42; Ag = 107.8682; Cd = 112.411; In = 114.82; Sn = 118.71;
Sb = 121.75; Te = 127.6; I = 126.90447; Xe = 131.29; Cs = 132.90543; Ba = 137.327; La = 138.9055;
Ce = 140.115; Pr = 140.90765; Nd = 144.24; Pm = 146.9151; Sm = 150.36; Eu = 151.965; Gd = 157.25;
Tb = 158.92534; Dy = 162.5; Ho = 164.93032; Er = 167.26; Tm = 168.93421; Yb = 173.04; Lu = 174.967;
Hf = 178.49; Ta = 180.9479; W = 183.85; Re = 186.207; Os = 190.2; Ir = 192.22; Pt = 195.08;
Au = 196.96654; Hg = 200.59; Tl = 204.3833; Pb = 207.2; Bi = 208.98037; Po = 208.9824; At = 209.9871;
Rn = 222.0176; Ac = 223.0197; Th = 226.0254; Pa = 227.0278; U = 232.0381; Np = 231.0359; Pu = 238.0289;
Am = 237.0482; Cm = 244.0642; Bk = 243.0614; Cf = 247.0703; Es = 247.0703; Fm = 251.0796; Md = 252.0829;
No = 257.0951; Lr = 258.0986; Rf = 259.1009; Db = 260.1053; Sg = 261.1087; Bh = 262.1138; Hs = 263.1182;
Mt = 262.1229;
Nn	= 1;	%Nn = Not named - for not named substances

if substance(1) ~= '"',         % substance chemical formula (not a name)

substance_index = 0;
bracket_index = 0;
last_bracket = [];
substance_char(1).second = '';
substance_char(1).first = '';
decimalFlag = 0;
while char_index < length(substance),
   char_index = char_index + 1;
   if (substance(char_index)>= 'A') & (substance(char_index)<= 'Z'),         
      substance_index = substance_index+1;												
      factor(substance_index) = 1;
      substance_mass(substance_index) = 0;
      substance_char(substance_index).first = substance(char_index);
   	number = 0;
   elseif (substance(char_index)>= 'a') & (substance(char_index)<= 'z'),		%small letter
      substance_char(substance_index).second = substance(char_index);
      number = 0;
   elseif substance(char_index) == '.' %start of decimal
       decimalFlag = 1;
   elseif (substance(char_index)>= '0') & (substance(char_index)<= '9'),		%number(could be decimal)
      number = 1;
      chartype(char_index,:) = 'nu';
      if decimalFlag == 0
          if chartype(char_index-1,:) == 'nu',
             factor(substance_index) = factor(substance_index)*10+str2num(substance(char_index));
          else%if 
             factor(substance_index) = str2num(substance(char_index));
          end%if
      else %deal with decimal
          factor(substance_index) = factor(substance_index) + str2num(substance(char_index))/(10^decimalFlag);  
          if char_index < length(substance)
                if (substance(char_index+1)>= '0') & (substance(char_index+1)<= '9')
                    decimalFlag = decimalFlag+1;
                else %end of decimal
                    decimalFlag = 0;
                end
          end
      end
   elseif (substance(char_index)== '+')|(substance(char_index)== '-'),
      factor(substance_index) = 1;      
   elseif substance(char_index)== '(',
      substance_index = substance_index+1;
      bracket_level = 1;
      bracket_char_index = 0;
      while bracket_level > 0,
         char_index = char_index+1;
         bracket_char_index = bracket_char_index+1;
         if 		substance(char_index) == '(', bracket_level = bracket_level + 1;
         elseif	substance(char_index) == ')', bracket_level = bracket_level - 1;
         end%if
         if bracket_level > 0,
            bracket_substance(bracket_char_index) = substance(char_index);
         end%if
      end%while
      substance_char(substance_index).first = '('; 
      substance_mass(substance_index) = MolMass(bracket_substance);
      number = 0;
   end%if
end%while char_index

for ei = 1:substance_index,
   if substance_char(ei).first ~= '(',
      if isempty(substance_char(ei).second),
         substance_mass(ei) = eval(substance_char(ei).first);
   	else
      	substance_mass(ei) = eval([substance_char(ei).first;substance_char(ei).second]);
      end% if
   end%if 
end%for
MM = substance_mass*factor';

elseif any(substance=='('),         % substance is not a chemical formula but a name 
  for i = 2:length(substance),
      if substance(i) == '(',
          j = i+1; MM_str = '';
          while substance(j) ~= ')' & j <= length(substance),
             MM_str = [MM_str substance(j)];
             j = j+1;
          end % while
      MM = str2num(MM_str);
      end % if
  end % for
else
  MM = 1;
end % if substance(1) ~= '"',         % substance is not a chemical formula but a name 

end % if substance(1) == '[',		  % vector of substances

%end MolMass.m