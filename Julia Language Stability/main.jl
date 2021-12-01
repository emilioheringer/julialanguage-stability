#Desmarcar caso nao tenha a biblioteca PrettyTables, ControlSystems e DPS instalada
#include("configuracoes/configuracoes.jl");

#Bibliotecas Utilizadas
using SparseArrays
using LinearAlgebra
using PrettyTables
using DSP
using ControlSystems;
using Plots;
#using PyCall;
#@pyimport(numpy);
#@pyimport(slycot);
#@pyimport(matplotlib);
    

#Inclusao das bibliotecas criadas
include("inclusoes.jl");
##adicionar os dados do sistema
include("dados_sistema.jl");
##Fluxo de Potência
include("fluxo_potencia.jl");
##Condicoes iniciais
include("condicoesiniciais.jl");
PD[2]=PD[2]-PD[2]*2;;
Yl = Complex.(zeros(nbus));
Sl=Complex.(zeros(nbus));
for i=1:nbus
    Sl[i] = PD[i] + im*QD[i];
        Yl[i]=conj(Sl[i])/V[i]^2;
end
Yl=Diagonal(Yl);
Y_bus=Y_bus+Yl;
#PD[2]=PD[2]*1.2;;
##Linearizacao
##Retorna a matriz de estados e entrada do sistema
A,B,C,D=decomposicao(linearizacao(Y_bus));
#realimentacao(A,B,C,D);

