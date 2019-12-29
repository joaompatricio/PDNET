c     $Id: dnsubs.f,v 1.1.1.1 2002/09/25 18:11:29 ggv Exp $

      subroutine dnbfs(nodes, narcs, laymax, dnnode, dnsrce, dnsink, 
     +                 dnnop2,dnlist, dndist, dnfapt, dnptrf, dnbapt,
     +                 dnptrb, dnfadj, dnbadj, dncap , dnflow,
     +                 dnbtof)
c     
c     dimension these to at least max number of nodes + 2:     
c            
      integer  nodes, narcs

      integer  dnlist(1:nodes+2), dndist(1:nodes+2), dnfapt(1:nodes+2), 
     +         dnptrf(1:nodes+2), dnbapt(1:nodes+2), dnptrb(1:nodes+2)  
c
c     dimension these to at least max number of arcs + 1:
c                   
      integer  dnfadj(1:narcs+1), dnbadj(1:narcs+1), dncap (1:narcs+1), 
     +         dnflow(1:narcs+1), dnbtof(1:narcs+1)        
c     
      integer  dnnode, dnsrce, dnsink, dnnop2, i2    , k2    ,
     +         j2    , qhead2, qtail2, laymax, k2dst , k2dp1 ,
     +         maxd2 , ibeg4 , iend4 , i4    , jjdn  , layer2,
     +         layp12   
c     
c*******************************************************************
c     
c     calling conditions:                                                   
c     internal call from dnsub only.                                  
c     input arrays:                                                         
c     dnfapt, dnbapt, dnfadj, dnbadj, dnbtof, dncap, dnflow           
c     scratch arrays:                                                       
c     dnlist                                                          
c     output:                                                               
c     in calling sequence (integer*2):                                    
c     laymax: 0 (layered net complete with sink in last layer.)       
c     1 (last layer is empty; current flow in original        
c     network is maximum.)                                 
c     in common (arrays):                                                 
c     dnptrf, dnptrb, dndist                                          
c***********************************************************************
c     
c     this subroutine constructs the layered network of the residual      
c     network w.r.t. the current flow.  the residual network is           
c     inferred, and the layered network is only recorded by temporary     
c     markers on the original graph. the subroutine implements a          
c     breadth-first search on the residual network, i.e.  nodes are       
c     processed as in a queue.  such a search partitions the nodes into   
c     a set of 'layers', according to their cardnality distance from      
c     the source. the dndist(.) array is used to record these distances,  
c     and thus defines the layers (note: dndist(isourc)=0).               
c     
c     in selecting arcs of the residual network that are directed from    
c     one layer to the next, we shall mark the arcs in the original       
c     network as 'open' or 'closed'.  a 'closed' arc (v,w) in the         
c     layered network is either an arc i=(v,w) that is saturated in the   
c     original network, or an arc i=(w,v) of the original network with    
c     zero flow. in the former case, the absence of this arc from the     
c     layered network is marked by negating dnfadj(i), and in the latter  
c     case, by negating dnbadj(i).  an 'open' arc i=(v,w) for the layered 
c     network either has positive residual capacity, or arc (w,v) has     
c     positive flow, in the original network. if a node v has no open     
c     arcs in the layered network that are in the forward adjacency list  
c     of v in the original network, then we set dnptrf(v)=0; and if it has
c     no open arcs that are in the backward adjacency list of v in the    
c     original network, we set dnptrb(v)=0.                               
c     
c     in this bfs search we only use open arcs that lead us to new        
c     nodes, or to nodes that have already been placed in the next        
c     layer. initially, the queue contains only the source.  as the       
c     vertices are popped from the queue and scanned, new nodes are       
c     injected into the queue. eventually, either the sink is reached,    
c     or some vertices are not reachable with open arcs from the          
c     current layer. in the former case the layered network for the       
c     current stage is complete, and thus a flow augmentation is          
c     possible.  in the latter case, the current flow on the original     
c     network is maximum, and the run terminates.                         
c     
c     the queue is maintained in the list dnlist(.) with two pointers,    
c     as shown:                                                           
c                   -------------------------------------------         
c    array              |    nodes to be scanned    |                   
c    dnlist(.):  ...    |                           |                   
c                   -------------------------------------------         
c                      ^                           ^                    
c                      |                           |                    
c                    qhead2                     qtail2                  
c*******************************************************************    
c                                                                       
c initialize:
c
      layer2 = 0
      do 10 i2=1,dnnode                                       
         dnptrf(i2)=0                                    
         dnptrb(i2)=0                                    
         dndist(i2)=dnnop2
 10   continue
      dndist(dnsrce)=0 
c                                       
c     put source into queue:                                                
c
      qhead2=0                                                
      qtail2=1                                                
      dnlist(1)=dnsrce                                        
      maxd2=dnnode-1                                          
c     
c     scan each node in queue:                                              
c     
c---------------
c     
 20   if(qhead2.eq.qtail2) goto 90                                    
      qhead2=qhead2+1                                         
c     pop node in front of queue:                                           
      k2=dnlist(qhead2)                                       
      k2dst=dndist(k2)                                        
      if(k2dst.ge.maxd2) goto 20                                       
      k2dp1=k2dst+1                                           
c     
c     scan node k2, i.e. search over fwd adjacency of k2 and for arcs      
c     (k2,j2) such that j2 is unscanned and j2 is not in the queue and     
c     arc (k2,j2) has positive residual capacity in original network.      
c     'j2 unscanned' is checked by the condition: 'dndist(j2)>=dndist(k2)'.
c     'j2 not in queue' is checked by the condition                        
c     'dndist(j2)=dndist(k2)+1 for      a scanned j2'.                     
c     
      ibeg4=dnfapt(k2)                                        
      iend4=dnfapt(k2+1)-1                                    
      if(ibeg4.gt.iend4) goto 60                              
      do 50 i4=ibeg4,iend4                                   
         j2=dnfadj(i4)                                   
         if(dndist(j2).le.k2dst) goto 40                  
         if(dncap(i4).le.dnflow(i4)) goto 40              
c     (implicitly) mark arc (k2,j2) as open in layered network; also,       
c     save its index in dnptrf(.):                                          
         dnptrf(k2)=i4                                   
         if(j2.ne.dnsink) goto 30                         
c     j2=sink:                                                              
         maxd2=k2dp1                                     
         dndist(j2)=maxd2                                
         goto 50                                         
c     append node j2 to the queue, if not already there:                    
 30      if(dndist(j2).eq.k2dp1) goto 50                  
         dndist(j2)=k2dp1                                
         qtail2=qtail2+1                                 
         dnlist(qtail2)=j2                               
         goto 50                                         
c     (explicitly) mark arc (k2,j2) as closed in layered network:           
 40      dnfadj(i4)=-dnfadj(i4)                          
 50   continue                                                
c     
c---------------
c     
c     if sink was reached, don't need to scan backward arcs      into its la
 60   if(dndist(dnsink).ne.dnnop2) goto 20                           
c     
c     else,      continue scan over forward arcsout of k2 in the residual   
c     network (which are backward arcs of the original network with         
c     positive flow). this code segment is analogous to the      one above: 
c     
      ibeg4=dnbapt(k2)                                                  
      iend4=dnbapt(k2+1)-1                                              
      if(ibeg4.gt.iend4) goto 20                                         
      do 80 i4=ibeg4,iend4                                            
         j2=dnbadj(i4)                                                   
         if(dndist(j2).le.k2dst) goto 70                                 
         jjdn=dnbtof(i4)                                                 
         if(dnflow(jjdn).le.0) goto 70                                   
         dnptrb(k2)=i4                                                   
         if(dndist(j2).eq.k2dp1) goto 80                                 
         dndist(j2)=k2dp1                                                
         qtail2=qtail2+1                                                 
         dnlist(qtail2)=j2                                               
         goto 80                                                        
 70      dnbadj(i4)=-dnbadj(i4)                                      
 80   continue                                                     
      goto 20                                                           
c     
c---------------
c     
 90      layer2=layer2+1                                               
         layp12=layer2+1                                                   
c     here,      all nodes in queue havebeen processed.                     
c     if the sink was reached, the layered network is complete:             
         laymax=0                                                          
         if(dndist(dnsink).ne.dnnop2)return                                
c     else,      the current flow is maximum:                               
         laymax=1                                                          
         return                                                            
c     
         end                                                               
      subroutine dnclea(dnarc,dnfadj,dnbadj)                            
c     
c     dimension these to at      least max number of arcs + 1:              
      integer dnarc,i4                                                  
      integer dnfadj(1:dnarc),dnbadj(1:dnarc) 
c     
c     
c*******************************************************************    
c calling conditions:                                                   
c      internal call from dnsub only.                                   
c input      (in common):                                               
c   scalars: dnarc                                                      
c   arrays:  dnfadj, dnbadj                                             
c output:                                                               
c   arrays (in common):      dnfadj,dnbadj                              
c********************************************************************   
c                                                                       
c clear      markings in arrays dnfadj() and dnbadj().:                  
      do 10 i4=1,dnarc                                                 
         if(dnbadj(i4).lt.0)dnbadj(i4)=-dnbadj(i4) 
         if(dnfadj(i4).lt.0)dnfadj(i4)=-dnfadj(i4)
 10   continue                                                       
      return                                                            
      end                                                               
      subroutine dncut(nodes,narcs,nncut,nacut,dnnode,dnnop2,dndist,
     +                 dnfapt,dnfadj,dncap,dnflow)   
c     
      integer nodes,narcs,nncut,nacut,dnnode,dnnop2,i2,ii2,ibeg4,iend4,
     +        i4,jjdn      
c     
c     
c     dimension these to at least max number of nodes + 2:                  
      integer  dndist(1:nodes+2),dnfapt(1:nodes+2)
c     
c     dimension these to at least max number of arcs + 1:                   
      integer dnfadj(1:narcs+1),dncap(1:narcs+1),dnflow(1:narcs+1)
c***********************************************************************
c     calling conditions:                                                   
c      only after a successfulreturn from a call to dnsub.              
c input      arrays:                                                    
c      dnfapt,dnfadj, dncap,dnflow, dndist                              
c output arrays:                                                        
c      dnfapt,dnfadj                                                    
c passed in calling sequence:                                           
c      nncut:number of nodes on source side of final cut (integer*4)    
c      nacut:number of arcs in the final cut (integer*4)                
c***********************************************************************
c                                                                       
c negate dnfapt(.) for nodes on      sink side of cut:                  
c
      nncut=0                                                           
      do 20 i2=1,dnnode                                                
         if(dndist(i2).ne.dnnop2) goto 10                        
         dnfapt(i2)=-iabs(dnfapt(i2))                           
         goto 20                                               
 10      nncut=nncut+1                                                 
 20   continue                                                      
c
c     negate dnfadj(.) for those saturated arcs in min cut:                 
c
      nacut=0                                                           
      do 40 i2=1,dnnode                                                
         if(dnfapt(i2).lt.0) goto 40
c                                       
c     node i2 is on      source side:                                       
c
         ibeg4=iabs(dnfapt(i2))                               
         iend4=iabs(dnfapt(i2+1))-1                          
         do 30 i4=ibeg4,iend4                               
            jjdn=dnfadj(i4)                                
            if(dnfapt(jjdn).gt.0) goto 30
c                                   
c     node dnfadj(i4) on sink side:                                         
c
            if(dnflow(i4).ne.dncap(i4)) goto 30           
            ii2=dnfadj(i4)                               
            dnfadj(i4)=-iabs(ii2)                       
            nacut=nacut+1                              
 30      continue                                                     
 40   continue                                                     
      return                                                            
      end                                                               
      subroutine dndfs(nodes,narcs,dnlfva,dnlaug,dnibig,dnsrce,dnsink,
     +     dnlist,dnfapt,dnptrf,dnbapt,dnflab,dnptrb,dnfadj,dnbadj,
     +     dncap,dnflow,dnbtof)
c     
c                                                                       
      integer  dnlfva, dnlaug, dnibig, dnsrce, dnsink, k2    ,
     +         j2    , kp2   , i4    , ipt4  , ifl4  , jjdn  ,
     +         nodes , narcs
c                                                                       
c dimension these to at least max number of nodes + 2:
c                  
      integer  dnlist(1:nodes+2), dnfapt(1:nodes+2), dnptrf(1:nodes+2), 
     +         dnbapt(1:nodes+2), dnflab(1:nodes+2), dnptrb(1:nodes+2) 
c
c dimension these to at least max number of arcs + 1:                   
c
      integer  dnfadj(1:narcs+1), dnbadj(1:narcs+1), dncap(1:narcs+1), 
     +         dnflow(1:narcs+1), dnbtof(1:narcs+1)          
c                                                                       
c***********************************************************************
c                                                                       
c  this      subroutine finds a maximal flowin the layered network. it  
c  uses      a depth-first search starting from the source and records  
c  the dfs tree      by the predecessor array dnlist(.) and flow label  
c  array dnflab(.). layered network arcs that are forward in the        
c  original network are      searched beforebackward arcs. the          
c  forward and backward      adjacencies of a node k2 are scanned inlifo
c  order, starting from      the last arc that was visited during a     
c  previous scan of k2.      thus, with the aid of the two pointer lists
c  dnptrf(.) and dnptrb(.), the      adjacency listsof the original     
c  network are scanned only once by this subroutine. dnptrf(k2)=0       
c  implies that      node k2has no outgoing unscanned arcs; dnptrb(k2)=0
c  implies that      node k2has no incoming unscanned arcs; a 'closed'  
c  node      is detected by the condition dnptrf(k2)=dnptrb(k2)=0.      
c                                                                       
c***********************************************************************
c                                                                       
c calling conditions:                                                   
c      internal call from dnsub only.                                   
c input      arrays:                                                    
c      dnptrb,dnptr, dnbtof, dncap, dnflow, dnbadj, dnfadj              
c output arrays:                                                        
c      dnptrf,dnptrb, dnlist, dnflab.                                   
c***********************************************************************
c                                                                       
c      subroutines called:dnpush                                        
c                                                                       
c***********************************************************************
c                                                                       
      k2=dnsrce                                                         
      dnlfva=0                                                         
      dnlaug=0                                                          
      dnflab(k2)=dnibig                                                 
c                                                                       
c-----------                                                            
c                                                                       
c scan node k2:                                                         
c                                                                       
 10   if(dnptrf(k2).eq.0) goto 50                                    
c find an open arc from      k2 to some nodej2: scan the                
c forward adjacency list of k2,      starting from  dnptrf(k2) and      
c proceeding backward on dnfadj() toward dnfapt(k2):                    
      i4=dnptrf(k2)                                                     
      ipt4=dnfapt(k2)                                                   
 20   j2=dnfadj(i4)                                                  
      if(j2.ge.0) goto 40                                              
      i4=i4-1                                                    
      if(i4.ge.ipt4) goto 20 
c     no open arcs out of k2; mark k2 as closed:                            
      dnptrf(k2)=0                                                      
      goto 50                                                          
c     an open arc found; move pointer dnptrf() to this arc;                 
c     extend dfs tree to node j2:                                           
 40   dnptrf(k2)=i4                                                
      dnflab(j2)=dnflab(k2)                                             
      if(dncap(i4)-dnflow(i4).lt.dnflab(j2))                            
     *     dnflab(j2)=dncap(i4)-dnflow(i4)                           
      dnlist(j2)=k2                                                     
      k2=j2                                                             
c     if node k2 (formerly j2) is the sink,      augmentflow; get new k2:   
      if(k2.eq.dnsink)                                                  
     *     call dnpush(nodes,narcs,k2,dnlfva,dnlaug,dnsrce,dnsink,
     +                 dnlist,dnptrf,dnflab,dnptrb,dnfadj,dnbadj,dncap,
     +                 dnflow,dnbtof)               
      goto 10                                                           
c     
c-----------
c     
 50   if(dnptrb(k2).eq.0) goto 80                                  
c     scan backward      arcs into k2 tofind an open arc:                   
      i4=dnptrb(k2)                                                     
      ipt4=dnbapt(k2)                                                   
 60   j2=dnbadj(i4)                                                
      if(j2.ge.0)goto 70                                              
      i4=i4-1                                                          
      if(i4.ge.ipt4)goto 60                                           
c     no open arc into k2:                                                  
      dnptrb(k2)=0                                                      
      goto 80                                                          
c     an open arc found; extend dfs      tree tonode j2:                    
 70   dnptrb(k2)=i4                                                
      dnflab(j2)=dnflab(k2)                                             
      jjdn=dnbtof(i4)                                                   
      ifl4=dnflow(jjdn)                                                 
      if(ifl4.lt.dnflab(j2))dnflab(j2)=ifl4                             
      dnlist(j2)=-k2                                                    
      k2=j2                                                             
      goto 10                                                           
c     
c-----------                                                            
c                                                                       
c k2 is      'closed' node; if it isthe source, then we have a maximal  
c flow in layered network:                                              
 80   if(k2.eq.dnsrce)return                                        
c else,      back upone node from k2 in dfs tree:                       
      kp2=dnlist(k2)                                                    
      k2=kp2                                                            
      if(k2.lt.0)k2=-k2                                                 
      if(kp2.ge.0) goto 90                                              
      i4=dnptrb(k2)                                                     
      dnbadj(i4)=-dnbadj(i4)                                            
      goto 10                                                           
 90   i4=dnptrf(k2)                                                 
      dnfadj(i4)=-dnfadj(i4)                                            
      goto 10                                                           
c                                                                       
c-----------                                                            
      end                                                               
      subroutine dnfwd(nodes,narcs,iretn,dnlist,dnfapt,dnfadj,          
     *   dnfrom,dnto,dncap,dngcap,dnmap)                                
                                                                        
c 01/14/92 - hack by geraldo veiga                                      
c save a edges -> forward adjancency map in dnmap                       
c                                                                       
      integer nodes,narcs,iretn,i2,i4,itail2,itpi4,itpsv4,itput4,jjdn,  
     *   nodp1 
                                                         
c dimension these to at least max number of nodes + 2:                  
      integer  dnlist(1:nodes+2),dnfapt(1:nodes+2)
c                                                                       
c dimension these to at least max number of arcs + 1:                   
      integer dnfrom(1:narcs+2),dnfadj(1:narcs+2),dnto(1:narcs+2),
     +        dncap(1:narcs+2),dngcap(1:narcs+2),dnmap(1:narcs+2)   
c                                                                       
c                                                                       
c***********************************************************************
c                                                                       
c calling conditions:                                                   
c      user callable; before calling dnsub.                             
c input:                                                                
c   in calling sequence      (all integer*4):                           
c      nodes=number of nodes, including source and sink                 
c      narcs=number of arcs.                                            
c                                                                       
c   in common, provide three arc lists with arcs in any      order:     
c      dnfrom=the list of tails. (integer*2)                            
c      dnto=the list of heads. (integer*2)                              
c      dngcap=the list of capacities. (integer*4).                      
c                                                                       
c scratch array      (in common):                                       
c      dnlist                                                           
c output:                                                               
c   in calling sequence      (integer*4):                               
c      iretn=0 no errors;  =1 error detected, check input.              
c   in common, arrays:                                                  
c      dnfapt,dnfadj,dncap (i.e. the fwd adjacency input data structure 
c       that is required as input to subroutine dnsub).                 
c                                                                       
c***********************************************************************
c                                                                       
c  using the three arc lists dnfrom(.),      dnto(.), dngcap(.), this   
c  subroutine constructs the forward adjacency arrays dnfapt(.),        
c  dnfadj(.), and dncap(.), as the input required by subroutine      dns
c                                                                       
c***********************************************************************
c                                                                       
c initialize:                                                           
      iretn=0                                                           
      nodp1=nodes+1                                                     
      do 10 i2=1,nodp1                                                 
         dnfapt(i2)=0                                                      
         dnlist(i2)=0
 10   continue
c temporarily store in dnfapt(.) number      of arcsout of each node:   
      do 20 i4=1,narcs                                                 
         jjdn=dnfrom(i4)                                                 
         dnfapt(jjdn)=dnfapt(jjdn)+1
 20   continue
c construct dnfapt(.):                                                  
      itpi4=dnfapt(1)                                                   
      dnfapt(1)=1                                                       
      do 30 i2=1,nodes                                                 
         itpsv4=itpi4+dnfapt(i2)                                           
         itpi4=dnfapt(i2+1)                                                
         dnfapt(i2+1)=itpsv4
 30   continue
c construct dnfadj(.) and dncap(.):                                     
      do 40 i4=1,narcs                                                 
         itail2=dnfrom(i4)                                               
         itput4=dnfapt(itail2)+dnlist(itail2)                            
         if(itput4.le.0) goto 50                                          
         dnfadj(itput4)=dnto(i4)                                         
         dnlist(itail2)=dnlist(itail2)+1                                 
         dncap(itput4)=dngcap(i4)                                        
         dnmap(i4)=itput4                                                
 40   continue                                                          
      return                                                            
 50   iretn=1                     
      return                                                            
      end                                                               
      subroutine dnpush(nodes,narcs,k2,dnlfva,dnlaug,dnsrce,dnsink,
     +           dnlist,dnptrf,dnflab,dnptrb,dnfadj,dnbadj,dncap,dnflow,
     +           dnbtof)               
c                                                                       
      integer nodes,narcs,dnlfva,dnlaug,dnsrce,dnsink,k2,j2,ksat2,
     +        incre4,i4,ii4     
c                                                                       
c dimension these to at least max number of nodes + 2:                  
      integer dnlist(1:nodes+2),dnptrf(1:nodes+2),dnflab(1:nodes+2),
     +        dnptrb(1:nodes+2)                   
c                                                                       
c dimension these to at least max number of arcs + 1:                   
      integer dnfadj(1:narcs+1),dnbadj(1:narcs+1),dncap(1:narcs+1),
     +        dnflow(1:narcs+1),dnbtof(1:narcs+1)          

c***********************************************************************
c                                                                       
c augment flow along flow augmenting path defined subroutine dndfs, i.e.
c using the predecessorarray dnlist(.), start from the sink and         
c traverse to the source arc flows and flow labels by an amount      equ
c to dnflab(dnsink). mark saturated forward (from s to t) arcs,      and
c backward arcs      having zero flow, as closed.                       
c                                                                       
c***********************************************************************
c                                                                       
c calling conditions:                                                   
c      internal call from dndfs only.                                   
c input      arrays:                                                    
c      dnflab,dnlist,dnptrb,dnptrf,dnbtof,dncap,dnflow,dnbadj,dnfadj    
c output arrays:                                                        
c      dnflab,dnflow,dnbadj,dnfadj                                      
c***********************************************************************
c                                                                       
      ksat2=0                                                           
      j2=dnsink                                                         
      incre4=dnflab(j2)                                                 
      dnlfva=dnlfva+incre4                                              
      dnlaug=dnlaug+1                                                   
c                                                                       
 10   k2=j2                                                          
      j2=dnlist(k2)                                                     
      dnflab(k2)=dnflab(k2)-incre4                                      
      if(j2.gt.0) goto 20  
c                                                                       
c decrease flow      on backward arc(k2,j2):                            
      j2=-j2                                                            
      i4=dnptrb(j2)                                                     
      ii4=dnbtof(i4)                                                    
      dnflow(ii4)=dnflow(ii4)-incre4                                    
      if(dnflow(ii4).ne.0) goto 10  
c flow is zero;      mark arc as 'closed' inlayered network:            
      dnbadj(i4)=-dnbadj(i4)                                            
      ksat2=j2                                                          
      goto 10                                                           
c                                                                       
c increase flow      on forward arc (j2,k2):                            
20      i4=dnptrf(j2)                                                  
      dnflow(i4)=dnflow(i4)+incre4                                      
      if(dncap(i4).ne.dnflow(i4)) goto 30   
c arc is now saturated;      mark isas 'closed' in layered net          
      dnfadj(i4)=-dnfadj(i4)                                            
      ksat2=j2                                                          
 30   if(j2.ne.dnsrce) goto 10                                        
c                                                                       
c return to resume search from k2 (node      closestto source of        
c the arc closed last)                                                  
      if(ksat2.eq.0) return                                               
      k2=ksat2                                                          
      return                                                            
      end                                                               
      subroutine dnsub(nodes,narcs,isourc,isink,maxflo,numaug,numstg,   
     *   nncut,nacut,iretn,dnlist,dndist,dnfapt,dnptrf,dnbapt,    
     *   dnflab,dnptrb,dnfadj,dnbadj,dncap,dnflow,   
     *   dnbtof)                                                        
c                                                                       
c************start dnsub standard user xface declaration block********* 
c                                                                       
c      user must include this block inhis/her calling program.          
c                                                                       
c   these declarations account for the entire array storage used by     
c   the dnsub subroutines, including what is necessary to storethe      
c   input data. it amounts to 6*nodes +5*arcs words, but it may be      
c   stated more precisely as follows:                                   
c                                                                       
c   1) for solving problems with up to 32,765 nodes and 2,147,483,647   ZZ
c   arcs, and with flow values of up to2,147,483,647, it suffices       
c   to declare arrays as shown under 'std arrangement' in the           
c   table below.  this arrangement is the one used in the release       
c   version of the dnsub subroutines.                                   
c                                                                       
c   2) in order to handle problems withmore than 32,765 nodes, all      
c   integer*2 declarations below must be changed to integer*4 and       
c   the      dnsub subroutines and user calling programs must be compile
c   and      linked.(see alternate 1 below.)                            
c                                                                       
c   3) for solving only      problems with at most 32,765 nodes, 32,767 
c   arcs, and with flows not exceeding 32,767, all integer*4            
c   declarations below may be changed to integer*2, and      the dnsub  
c   subroutines       and user calling programs mustbe compiled and     
c   linked. (see alternate 2 below.)  it is also necessary to           
c   activate the statement dnibig=32767      in subroutine dnsub, and to
c   include 'implicit integer*2 (i-n)' statements in each program       
c   unit. this  arrangement, which usesthe least amount of memory       
c   and      is the fastest when executed on16-bit processors, is only  
c   useful for a set of      very special applications.  itsgeneral     
c   use      is not recommended.  other alternates are also possible.   
c                                                                       
c  -----------------------------------------------------------------    
c                                                                       
c        std arrangement    alternate 1      alternate 2                
c             (not recommended)                                         
c        ---------------  ---------------  ----------------             
c  array   length type eqv'd type  type eqv'd type  type eqv'd type     
c  ------  -----    ----  -------  ---- ----- ----  ---- ----- -----    
c  dnfapt    n      4               4                2                  
c  dnfadj    a      2               4                2                  
c  dncap     a      4               4                2                  
c  dnbapt    n      4               4                2                  
c  dnbadj    a      2 dnfrom  2     4 dnfrom  4      2 dnfrom   2       
c  dnflow    a      4 dnto    2     4 dnto    4      2 dnto     2       
c  dnbtof    a      4 dngcap  4     4 dngcap  4      2 dngcap   2       
c  dnptrf    n      4               4                2                  
c  dnptrb    n      4               4                2                  
c  dnlist    n      2               4                2                  
c  dnflab    n      4 dndist  2     4 dndist  4      2 dndist   2       
c  -----------------------------------------------------------------    
c                                                                       
c  total bytes:     22*n + 16*a     24*n + 20*a       12*n + 10*a       
c  (n=nodes, a=arcs)                                                    
c  -----------------------------------------------------------------    
c                                                                       
c subroutine arguments                                                  
      integer nodes,narcs,isourc,isink,maxflo,numaug,numstg,nncut,nacut,
     *   iretn                                                          
c dimension these to at least max number of nodes + 2:                  
      integer dnlist(1:nodes+2),dndist(1:nodes+2),dnfapt(1:nodes+2),
     +        dnptrf(1:nodes+2),dnbapt(1:nodes+2),dnflab(1:nodes+2),
     +        dnptrb(1:nodes+2) 
c                                                                       
c dimension these to at least max number of arcs + 1:                   
      integer dnfadj(1:narcs+1), dnbadj(1:narcs+1),dncap(1:narcs+1),
     +   dnflow(1:narcs+1),dnbtof(1:narcs+1)
                                                                        
c all arrays are passed to dnsub as arguments.  the original            
c equivalence statement can by simulated by passing identical pointers  
c in the calling sequence for each member of the equivalence group.     
c                                                                       
c      equivalence (dnflow(1),dnto(1)),(dnflab(1),dndist(1)),           
c     *              (dnfrom(1),dnbadj(1)),(dngcap(1),dnbtof(1))         
c                                                                       
c                                                                       
c declarations for original dn00 common                                 
      integer dnarc,dnfva,dnaug,dnlfva,dnlaug,dnstge,                   
     *   dnibig,dnnode,dnsrce,dnsink,dnnop2                             
c                                                                       
c local variables to this subroutine                                    
      integer  nodp12,i2,ihead2,laymax,jjdn,ihput4,i4,                  
     *   ibeg4,iend4,ihpi4,ihpsv4                                       
c                                                                       
                                                                        
c                                                                       
c************end dnsub standard user xface declaration block*********** 
c***********************************************************************
c                                                                       
c calling conditions:                                                   
c       user callable with input data as follows.                       
c input:                                                                
c   scalars (in calling sequence; all integer*4):                       
c       nodes:  number of nodes (including source and sink)             
c       narcs:  number of arcs                                          
c       isourc: node number for source                                  
c       isink:  node number for sink                                    
c                                                                       
c   arrays (in common):                                                 
c       dnfapt: (nodes+1)-long integer*4 pointer array for forward adja-
c               cency lists (i.e. the forward adjacency list of a node i
c               is the set of arcs                                      
c                      {(i,dnfadj(j)) :j=dnfapt(i),...,dnfapt(i+1)-1 }  
c               note: must have dnfapt(nodes+1) = narcs+1.              
c                                                                       
c       dnfadj: narcs-long integer*2 array giving the list of nodes in  
c               the forward adjacency list dnfadj(j) for each node j,   
c               as described above.                                     
c                                                                       
c       dncap:  narcs-long integer*4 array giving the arc capacities, in
c               the order prescribed by dnfapt(.) and dnfadj(.). all arc
c               capacities must be given as positive integers.          
c                                                                       
c   note:       for transforming unordered arc lists to the above input 
c               data structure,use subroutine dnfwd before calling dnsub
c                                                                       
c output:                                                               
c                                                                       
c   scalars (in calling sequence):                                      
c       maxflo: value of maximum flow (integer*4).                      
c       numaug: number of flow augmentations (integer*4).               
c       numstg: number of stages (layered networks created; integer*4). 
c       nncut:  number of nodes on source side of final cut (integer*4) 
c       nacut:  number of saturated arcs in the final cut (integer*4)   
c       eltim:  execution time, in seconds (see note 3 below) (real*4). 
c       iretn:  nonzero if there are errors in input data.              
c                                                                       
c   arrays (in common:):                                                
c       dnfapt: the original dnfapt(.) with some of its elements        
c               negated to mark the nodes that are on the sink side of  
c               the minimum cut, i.e. dnfapt(i)<0 implies that node i is
c               on the sink side of cut, else, it is on the source side.
c                                                                       
c       dnfadj: the original dnfadj(.) with some of its elements negated
c               to mark those saturated arcs that are in the min cut    
c               found by the algorithm, i.e.  dnfadj(i)<0 implies that  
c               the i-th arc in forward adjaceny order is in the cut;   
c               the following code segment prints these arcs:           
c                                                                       
c                               do 2 k=1,nodes                          
c                                  ibeg=iabs(dnfapt(k))                 
c                                  iend=iabs(dnfapt(k+1))-1             
c                                  do 1 i=ibeg,iend                     
c                                    if(dnfadj(i).lt.0)print the arc    
c                       1          continue                             
c                       2       continue                                
c                                                                       
c               (see note 1 below).                                     
c                                                                       
c       dncap:  the original capacities, unaltered.                     
c                                                                       
c       dnflow: narcs-long integer*4 array which gives flows on arcs, in
c               forward adjacency order.                                
c                                                                       
c       dnbapt: (nodes+1)-long integer*4 pointer array for backward     
c               adjacency lists (i.e. the forward adjacency list of     
c               node i is the set of arcs                               
c                   {(dnbadj(j),j): j=dnbapt(i),...,dnbapt(i+1)-1 }.    
c               note: we must have dnbapt(nodes+1) = narcs+1.           
c                                                                       
c       dnbadj: narcs-long integer*2 array giving the backward adjacency
c               list of each node j, one after the other, as described  
c               above (also see note 2 below).                          
c                                                                       
c       dnbtof: narcs-long integer*4 array giving the backward adjacency
c               to forward adjacency mapping, i.e. the j-th arc         
c               in the backward adjacency order is the dnbtof(j)-th     
c               arc in the forward adjacency order.                     
c                                                                       
c***********************************************************************
c                                                                       
c internal scratch arrays (in common):                                  
c                                                                       
c       dnlist: (nodes+1)-long integer*2 array used as follows:         
c               1) as the queue in bfs search for constructing the      
c                  layered graph (subr. dnbfs).                         
c               2) as the "parent" array in dfs search of layered       
c                  graph (subr: dndfs, dnpush).                         
c               3) as scratch in constructing fwd adjacencies           
c                  (subr dnfwd), and bwd adjacencies (subr. dnsub).     
c                                                                       
c       dnflab: nodes-long integer*4 array used in dfs search of        
c               layered graph to store the flow label for each node     
c               (see array dndist() ).                                  
c                                                                       
c       dndist: nodes-long integer*2 array used in constructing the     
c               layered graph by bfs, storing the distance of each      
c               node from the source; dndist(isourc)=0 (equivalenced    
c               to the integer*4 array dnflab() ).                      
c                                                                       
c       dnptrf: nodes-long integer*4 status and pointer array, i.e.     
c               dnptrf(k)=0 if there are no open arcs leaving node k    
c                               in the representation of the layered    
c                               network.                                
c                        >0 the arc index, in 'fa' order, that          
c                               corresponds to the last arc in the      
c                               fwd adjacency of k which is scanned     
c                               backward.                               
c                                                                       
c       dnptrb: nodes-long integer*4 status and pointer array, i.e.     
c               dnptrb(k)=0 if there are no open arcs entering node k   
c                               in the representation of the layered    
c                               network.                                
c                        >0 the arc index (in 'ba' order) of the first  
c                               arc in the bwd adjacency of k to be     
c                               scanned in lifo order.                  
c                                                                       
c***********************************************************************
c                                                                       
c notes:1) dnfadj(j) is internally negated to mark the j-th arc (in     
c          forward adjacency order) as 'closed'. before returning       
c          to the calling pgm, these markings are discarded, and some   
c          elements of dnfadj(.) are negated to reflect the output      
c          specification described above.                               
c                                                                       
c       2) dnbadj(j) is internally negated to mark the j-th arc         
c          (in backward adjacency order) as 'closed'.                   
c                                                                       
c       3) the user must provide a timing subroutine that returns       
c          the current time in variable ' t ', in seconds.  such        
c          a subroutine is installation-dependent and is not given here.
c          if this cannot be done, provide the following subroutine and 
c          forget the execution timing information:                     
c                                                                       
c               subroutine getime(t)                                    
c               real t                                                  
c               t=0.0                                                   
c               return                                                  
c               end                                                     
c                                                                       
c***********************************************************************
c                                                                       
c       subroutines called:     dnbfs, dndfs, dncut, getime,dnclea      
c                                                                       
c***********************************************************************
c                                                                       
c   names collectively reserved by all dnsub subroutines:               
c                                                                       
c   dn00, dn01, dn02, dn03, dn04, dn05, dn06, dn08, dn09, dn10, dn11,   
c   dn12, dnarc, dnaug, dnbadj, dnbapt, dnbfs, dnbtof, dncap, dnclea,   
c   dncut, dndfs, dndist, dnfadj, dnpush, dnfapt, dnflab, dnflow,       
c   dnfrom, dnfva, dnfwd, dngcap, dnibig, dnlaug, dnlfva, dnlist,       
c   dnnode, dnnop2, dnout, dnptrb, dnptrf, dnsink, dnsrce, dnstge,      
c   dnsub, dneltm, dnto.                                                
c                                                                       
c********************************************************************** 
c initialization:                                                       
c                     

        iretn=0                                                         
        dnibig=2147483647                                               
c        dnibig=32767                                                    
        dnarc=narcs                                                     
        dnnode=nodes                                                    
        dnnop2=dnnode+2                                                 
        dnsrce=isourc                                                   
        dnsink=isink                                                    
        do 10 i4=1,dnarc                                               
           dnflow(i4)=0
 10     continue
c***********************************************************************
c                                                                       
c using the forward adjacency arrays dnfapt(.) and dnfadj(.) (in        
c common), create the backward adjacency arrays dnbapt(.), dnbadj(.)    
c and the backward to forward adjacency mapping dnbtof(.).  dnlist(.) is
c used here as a scratch array.                                         
c              
        nodp12=dnnode+1                                                 
        do 20 i2=1,nodp12                                              
           dnlist(i2)=0                                            
           dnbapt(i2)=0                                            
 20     continue
c temporarily store in dnbapt(.) number of arcs into each node:         
        do 40 i2=1,dnnode                                              
           ibeg4=dnfapt(i2)                                        
           iend4=dnfapt(i2+1)-1                                    
           do 30 i4=ibeg4,iend4                                   
              jjdn=dnfadj(i4)                                       
              dnbapt(jjdn)=dnbapt(jjdn)+1 
 30        continue
 40     continue                      
c construct dnbapt(.):                                                  
        ihpi4=dnbapt(1)                                                 
        dnbapt(1)=1                                                     
        do 50 i2=1,dnnode                                              
           ihpsv4=ihpi4+dnbapt(i2)                                 
           ihpi4=dnbapt(i2+1)                                      
           dnbapt(i2+1)=ihpsv4 
 50     continue
c construct dnbadj(.) and dnbtof(.):                                    
        do 70 i2=1,dnnode                                              
           ibeg4=dnfapt(i2)                                        
           iend4=dnfapt(i2+1)-1                                    
           do 60 i4=ibeg4,iend4                                   
              ihead2=dnfadj(i4)                               
              ihput4=dnbapt(ihead2)+dnlist(ihead2)            
              if(ihput4.le.0) goto 100
              dnbadj(ihput4)=i2                               
              dnbtof(ihput4)=i4                               
              dnlist(ihead2)=dnlist(ihead2)+1 
 60        continue
 70      continue   
c                                                                       
c***********************************************************************
c                                                                       
        dnfva=0                                                         
        dnaug=0                                                         
        dnstge=1                                                        
c                                                                       
c stage loop:                                                           
 80     continue                                                        
c   form new layered network w.r.t. current feasible flow.              
                call dnbfs(nodes,narcs,laymax,dnnode,dnsrce,dnsink,
     +                     dnnop2,dnlist,dndist,dnfapt,dnptrf,dnbapt,
     +                     dnptrb,dnfadj,dnbadj,dncap,dnflow,dnbtof)
                if(laymax.gt.0) goto 90                                 
c   compute a maximal flow in this layered net:                         
                call dndfs(nodes,narcs,dnlfva,dnlaug,dnibig,dnsrce,
     +                     dnsink,dnlist,dnfapt,dnptrf,dnbapt,dnflab,
     +                     dnptrb,dnfadj,dnbadj,dncap, dnflow,dnbtof)
c   update flow value:                                                  
                dnfva=dnfva+dnlfva                                      
                dnaug=dnaug+dnlaug                                      
                dnstge=dnstge+1                                         
c   erase markings that define layered network:                         
                call dnclea(dnarc,dnfadj,dnbadj)                        
        goto 80                                                         
c***********************************************************************
c                                                                       
c current flow is maximum:                                              
 90     maxflo=dnfva                                                    
        numaug=dnaug                                                    
        numstg=dnstge                                                   
        call dnclea(dnarc,dnfadj,dnbadj)                                
c                                                                       
c mark the min cut-set in dnfapt(.) and dnfadj(.):                      
        call dncut(nodes,narcs,nncut,nacut,dnnode,dnnop2,dndist,dnfapt,
     +             dnfadj,dncap,dnflow)
c                                                                       
        return                                                          
c                                                                       
 100    iretn=1                                                         
        return                                                          
        end                                                             
