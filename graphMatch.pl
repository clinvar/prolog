
% outputMatchNodeListsWithGraphs(NodeLists, Graphs)
% Performs pairwise comparisons of lists of nodes against graphs. 
% Input: node/2 and edge/3 assertions included in this file
% Output: tabulated (rows: graphs; columns: nodes and totals) in file NGmatches.txt
% Example: outputMatchNodeListsWithGraphs([[n3,n2,n1],[n4,n5,n6]], [g0,g1,g2,g3]).
%
% Author: Aleks Milosavljevic 9/1/2014, BRL, BCM
%
% Annotated by Xin Feng, 09/04/2014
%

outputMatchNodeListsWithGraphs(NodeLists, Graphs):-
  % open the file for writing
  % `append` means writing
open('NGmatches.txt', append, OutStream),
matchNodeListsWithGraphs(NodeLists, Graphs, OutStream).



 % Input: []: An empty list
 %         _: Any variable
 %
matchNodeListsWithGraphs( [], _, OutStream):-
 % If the input list is empty, output a new line
nl(OutStream),
 % and then close the output stream
close(OutStream).

 % Input: [NodeList|RestNodeLists] : a list divided into
 %         NodeList and RestNodeLists
 %  
matchNodeListsWithGraphs( [NodeList|RestNodeLists], Graphs, OutStream):-
  % Output the header 
writeNodeListHeader(NodeList, OutStream),
  % Match the first nodelist within the Graphs
matchSingleNodeListWithGraphs( NodeList, Graphs, OutStream),
  % Match the rest nodelists 
matchNodeListsWithGraphs( RestNodeLists, Graphs, OutStream).

 % If the input nodelist is empty
matchSingleNodeListWithGraphs( [], _, OutStream):-
  % then just output a new line
nl(OutStream).

 % If the input Graph-to-match is empty
matchSingleNodeListWithGraphs( _, [], OutStream):-
  % then just output a new line
nl(OutStream).

matchSingleNodeListWithGraphs( NodeList, [NextGraph|RestGraphs], OutStream):-
 % Otherwise, start the matching with the first graph (NextGraph)
matchSingleNodeListWithSingleGraph( NodeList, NextGraph, OutStream),
 % Repeat for the rest of the graph
matchSingleNodeListWithGraphs( NodeList, RestGraphs, OutStream).

 % How to we actually match a node list with a grpah?
matchSingleNodeListWithSingleGraph(NodeList, Graph, OutStream):-
  % The match is determined by the information of the nodelist
relInfo(0.9,0.1,[],NodeList,Graph,L,T),
  % Output the match result
writeNodeListValues(Graph, T, L,  Graph, OutStream).


% I/O utilities
writeNodeListHeader(NodeList, OutStream):-
write( OutStream, 'Graph '),
write( OutStream, 'TotalRelInfo '),
writeNodeList(NodeList,OutStream).

writeNodeList([],OutStream):-
nl(OutStream).

writeNodeList([Node|RestNodes],OutStream):-
write(  OutStream, Node ),
write(  OutStream, ' ' ),
writeNodeList(RestNodes,OutStream).

writeNodeListValues(Graph, T, L,  Graph, OutStream):-
write( OutStream, Graph),
write( OutStream, ' ' ),
write( OutStream, T),
write(  OutStream, ' '),
writeNodeRelInfo(L,OutStream).

writeNodeRelInfo([],OutStream):-
nl(OutStream).

writeNodeRelInfo([X|LX],OutStream):-
write(  OutStream, X),
write(  OutStream,' '),
writeNodeRelInfo(LX,OutStream).

% graph utility functions

nodeOf(X,G):-
node(X,G).

nodeOf(X,G):-
edge(X,_,G).

nodeOf(X,G):-
edge(_,X,G).

allNodesInGraph(G,Lnodes):-
setof(X,nodeOf(X,G),Lnodes).

% Not sure about how this works
notNodeOf(X,G):-
nodeOf(X,G),
!,
fail.

% Not sure about this
notNodeOf(_,_).

nodeNotInGraph([X|_],G):-
notNodeOf(X,G).

nodeNotInGraph([X|T],G):-
nodeOf(X,G),
nodeNotInGraph(T,G).

nodeNotInGraph([],_):-
!,
fail.

%relInfo: calculate information for given nodes LX relative to graph G 
%NOTE: calculate information of all subsets of LX, in order of LX
%example1 invocation: relInfo(0.9,0.1,[],[n1,n2,n3],g0,L,T).
%example2 invocation: relInfo(0.9,0.1,[],[n1,n2,n3],g1,L,T).
%example3 invocation: relInfo(0.9,0.1,[],[n1,n2,n3],g2,L,T).
%example4 invocation: relInfo(0.9,0.1,[],[n1,n2,n3],g3,L,T).
%example5 invocation: relInfo(0.9,0.1,[],[n4,n5,n6],g0,L,T).
%example6 invocation: relInfo(0.9,0.1,[],[n4,n5,n6],g1,L,T).
%example7 invocation: relInfo(0.9,0.1,[],[n4,n5,n6],g2,L,T).
%example8 invocation: relInfo(0.9,0.1,[],[n4,n5,n6],g3,L,T).

relInfo(_,_,_,[],_,[],0).

relInfo(Pdiff,Peq,LXcurrent,[X|TLXrest],G,[XrelInfo|TLrelInfo],TotRelInfoPlus):-
relInfoNext(Pdiff,Peq,LXcurrent,G,X,XrelInfo),
append(LXcurrent,[X],LXplus),
relInfo(Pdiff,Peq,LXplus,TLXrest,G,TLrelInfo,TotRelInfo),
TotRelInfoPlus is XrelInfo + TotRelInfo.

relInfoNext(Pdiff,Peq,LX,G,X,XrelInfo):-
condP(Pdiff,Peq,LX,G,LP),
extractRelInfoNext(LP,X,XrelInfo).

extractRelInfoNext([],_,XrelInfo):-
maxNodeInfo(XrelInfo).

extractRelInfoNext([(_,P)|_],_,XrefInfo):-
P =< 0,
maxNodeInfo(XrelInfo).

extractRelInfoNext([(X,P)|_],X,XrefInfo):-
P > 0,
log(P,MinusXrefInfoInNats),
log(2,LogOf2),
XrefInfo is - (MinusXrefInfoInNats / LogOf2).

extractRelInfoNext([(N,_)|TLP],X,XrefInfo):-
not(compare(=,N,X)),
extractRelInfoNext(TLP,X,XrefInfo).

maxNodeInfo(1000).

%condP: given nodes LX and graph G, assign probabilities Pdiff+Peq =1 to remaining nodes in G 

condP(Pdiff,Peq,_,_,[]):-
Ptot is Pdiff + Peq,
1.0 < Ptot, 
!. %error

condP(Pdiff,Peq,_,_,[]):-
Ptot is Pdiff + Peq,
1.0 > Ptot, 
!. %error

condP(_,_,LX,G,[]):-
nodeNotInGraph(LX,G),
!. %error

condP(_,Peq,[],G,LP):-
eqWeights(Peq,G,LWeq),
addUpWeights(LWeq,LWtot),
normalizeWeights(LWtot,LP).

condP(Pdiff,Peq,LX,G,LP):-
length(LX,LXsize),
LXsize > 0,
Wdiff is Pdiff/LXsize,
diffWeights(Wdiff,LX,G,LWdiff),
eqWeights(Peq,G,LWeq),
append(LWdiff,LWeq,LW),
addUpWeights(LW,LWtot),
normalizeWeights(LWtot,LN),
collectRemainingNodes(LX,LN,LR),
normalizeWeights(LR,LP).

collectRemainingNodes(_,[],[]).

collectRemainingNodes(LX,[(N,_)|TLN],TLR):-
member(N,LX),
collectRemainingNodes(LX,TLN,TLR).

collectRemainingNodes(LX,[(N,W)|TLN],[(N,W)|TLR]):-
not(member(N,LX)),
collectRemainingNodes(LX,TLN,TLR).

diffWeights(_,[],_,[]).

diffWeights(Wdiff,[X|T],G,L):-
divideWeight(X,G,Wdiff,LX),
diffWeights(Wdiff,T,G,LT),
append(LX,LT,L).

eqWeights(Peq,G,L):-
allNodesInGraph(G,Lnodes),
length(Lnodes,NodeCount),
Pnode is Peq / NodeCount,
assignWeights(Pnode,Lnodes,L).

%example1 invocation: eqWeights(0.1,g3,L).
   
assignWeights(_,[],[]).

assignWeights(Pnode,[N|R],[(N,Pnode)|T]):-
assignWeights(Pnode,R,T).

addUpWeights(LW,LP):-
sort(LW,LS),
addUpSortedWeights(LS,LP).

addUpSortedWeights([],[]).

addUpSortedWeights([(N,Wnext)|T],[(N,Wtotal)|R]):-
addUpSortedWeightsForCurrentNode( N, [(N,Wnext)|T], TR, Wtotal ),
addUpSortedWeights(TR,R).

addUpSortedWeightsForCurrentNode( N, [(N,Wnext)|T], TR, Wtotal ):-
addUpSortedWeightsForCurrentNode( N, T, TR, Wcurrent ),
Wtotal is Wnext+Wcurrent.

addUpSortedWeightsForCurrentNode( _, [], [], 0 ).

addUpSortedWeightsForCurrentNode( N, [(Nnext,Wnext)|T], [(Nnext,Wnext)|T], 0 ):-
N \= Nnext.                                                            

normalizeWeights([],[]).

normalizeWeights(LP,LP):-
totalWeights(LP,0).

normalizeWeights(LP,LN):-
totalWeights(LP,Wtot),
Wtot > 0,
CorrectionFactor is 1 / Wtot,
multiplyWeightsByCorrectionFactor(CorrectionFactor,LP,LN).

totalWeights([],0).

totalWeights([(_,W)|R],Wtot):-
totalWeights(R,Wcurrent),
Wtot is Wcurrent + W.

multiplyWeightsByCorrectionFactor(_,[],[]).

multiplyWeightsByCorrectionFactor(CorrectionFactor,[(N,W)|T],[(N,Wcorrected)|R]):-
Wcorrected is CorrectionFactor*W,
multiplyWeightsByCorrectionFactor(CorrectionFactor,T,R).

subthresholdWeight(W):- 
W < 0.01.

%example1 invocation: divideWeight(n1,g1,0.9,L).

divideWeight(_,_,W,[]):-
subthresholdWeight(W).

divideWeight(X,G,W,[(X,W)]):-
findall(Y,edge(X,Y,G),NeighborsOfX),
length(NeighborsOfX,0).

divideWeight(X,G,W,[(X,Wx)|T]):-
findall(Y,edge(X,Y,G),NeighborsOfX),
length(NeighborsOfX,NumberOfNeighbors),
NumberOfNeighbors > 0,
Wx is W/2,
WneighborOfX is Wx/NumberOfNeighbors,
distributeWeightAmongNeighbors(NeighborsOfX,G,WneighborOfX,T).

distributeWeightAmongNeighbors(_,_,WneighborOfX,[]):-
subthresholdWeight(WneighborOfX).

distributeWeightAmongNeighbors([],_,_,[]).

distributeWeightAmongNeighbors([X|TX],G,W,T):-
divideWeight(X,G,W,Tfirst),
distributeWeightAmongNeighbors(TX,G,W,Trest),
append(Tfirst,Trest,T).

% graph 0 - contains disconneceted nodes
node(n1,g0).
node(n2,g0).
node(n3,g0).
node(n4,g0).
node(n5,g0).
node(n6,g0).

% graph 1 - contains clique {n1,n2,n3}
node(n4,g1).
node(n5,g1).
node(n6,g1).
edge(n1,n2,g1).
edge(n2,n1,g1).
edge(n1,n3,g1).
edge(n3,n1,g1).
edge(n2,n3,g1).
edge(n3,n2,g1).

% graph 2 - contains clique {n4,n5,n6}
node(n1,g2).
node(n2,g2).
node(n3,g2).
edge(n4,n5,g2).
edge(n5,n4,g2).
edge(n4,n6,g2).
edge(n6,n4,g2).
edge(n5,n6,g2).
edge(n6,n5,g2).

% graph 3 - contains bidirectional path n1<->n2<->n3<->n4<->n5<->n6
node(n1,g3).
node(n2,g3).
node(n3,g3).
node(n4,g3).
node(n5,g3).
node(n6,g3).
edge(n1,n2,g3).
edge(n2,n1,g3).
edge(n2,n3,g3).
edge(n3,n2,g3).
edge(n3,n4,g3).
edge(n4,n3,g3).
edge(n4,n5,g3).
edge(n5,n4,g3).
edge(n5,n6,g3).
edge(n6,n5,g3).




