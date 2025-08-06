function matrix=read_matrix(filename,nz,nx,format,machineformat)
% read a matrix from a file,
% matrix=read_matrix(filename,nz,nx,format,machineformat)
% filename: the file name which is read
% nz: the number of row; nx: the number of coloumn
% format: 'float32'; 'int32'
% machineformat : the format of the file
%                     'n' or 'native' The byte ordering that your system uses (default)
%                     'b' or 'ieee-be' Big-endian ordering
%                     'l' or 'ieee-le' Little-endian ordering
%                     's' or 'ieee-be.l64' Big-endian ordering, 64-bit data type
%                     'a' or 'ieee-le.l64' Little-endian ordering, 64-bit data type

if nargin<5
    machineformat='native';
end
if nargin<4
    format='float32';
end
fip=fopen(filename,'r');
matrix=fread(fip,[nz nx],format,machineformat);
fclose(fip);
end