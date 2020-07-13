%% Intitializations
clear;
clc;
prompt = {'Enter the file name of the network to be analyzed'};
title1 = 'Input';
dims =  [1 50];
definput = {'ko04141'};
answer = inputdlg(prompt,title1,dims,definput);
filename1=strcat('D:\rishika\4_work_boolean_logic\KEGG_data\',answer,'.xml');
filename=char(filename1);
fid=fopen(filename);
filename
line1=fgets(fid);
%line1
nodes=zeros(1,2);
%edges=[];
i=1;
temp=0;
while ischar(line1) && temp==0
       x1=strfind(line1,'relation entry1');
       if isempty(x1)==0
           %line1
           A1=strsplit(line1);
           p=strsplit(A1{3},'"');
           nodes(i,1)=str2double(cellstr(p{2}));
           p=strsplit(A1{4},'"');
           nodes(i,2)=str2double(cellstr(p{2}));
           line1=fgets(fid);
           if isempty(strfind(line1,'subtype'))==1
               
               continue;
           end
           A2=strsplit(line1);
           p=strsplit(A2{3},'"');
           %disp(cellstr(p{2}));
           edges(i,1)=cellstr(p{2});
           i=i+1;
           %line1
           
       end
       line1=fgets(fid);
       
end
fclose(fid);
g=unique(nodes);

node_list=nodes;
edge_list=edges;

%% Finding node numbers
fid=fopen(filename);
line1=fgets(fid);
i=1;
b=1;
temp=0;
while ischar(line1) && temp==0
       %line1
       x1=strfind(line1,'entry id');
       x2=strfind(line1,'path');
       
        if isempty(x1)==0 && isempty(x2)==1
            %line1
            %temp=1;
            A1=strsplit(line1);
            p=strsplit(A1{3},'"');
            label_num(i,1)=str2double(cellstr(p{2}));
            A2=strsplit(line1,'"');
            label_name(i,1)=cellstr(A2{4});
            
            line1=fgets(fid);
            x5=strfind(line1,'link');
            if isempty(x5)==0
              line1=fgets(fid);
            end
            x3=strfind(line1,'graphics name');
            x4=strfind(line1,'color');
            if isempty(x3)==0 && isempty(x4)==0
            
               A2=strsplit(line1,'=');
               A3=strsplit(A2{2}, '"');
               
               graphics_name(i,1)=cellstr(A3{2});
               
            end
            if isempty(x4)==0 && isempty(x3)==1
               %disp('entered');
               graphics_name(i,1)=cellstr('no name');
               
            end
            i=i+1;
        end
       line1=fgets(fid);
       
end
%text=nodes;

% edges
% nodes
j=1;k=1;
non_bind_edge={};
non_bind=[];

m=length(nodes);
for i=1:length(nodes)
    %temp=0;
    if strcmp(edges(i,1),'indirect')==0 && strcmp(edges(i,1),'missing')==0 && strcmp(edges(i,1),'compound')==0%&& strcmp(edges(i,1),'dephosphorylation')==0
        non_bind_edge(j,1)=edges(i,1);
        non_bind(j,1)=nodes(i,1);
        non_bind(j,2)=nodes(i,2);
        j=j+1;
               
    end
    
end
fprintf('Number of nodes %d \n',length(unique(non_bind)));
fprintf('Number of edges %d \n',length(non_bind_edge));
m1=non_bind;
a = unique(non_bind);
out = [a,histc(non_bind(:),a)];
%fprintf('Maximum connectivity %d \n',out(:,2));
f=max(out(:,2));
fprintf('Maximum connectivity of node %d \n',f);
f=sum(out(:,2) == f);
fprintf('Number of nodes with maximum connectivity %d \n', f);
temp_edge=edges;
temp_node=nodes;
edges={};
nodes=[];
edges=non_bind_edge;
nodes=non_bind;
%edges
length(nodes)

%% Initialization
initial_node_list=unique(nodes);
m=length(initial_node_list);
k=1;
for i=1:m
    initial_node_list(i,2)=k;
    k=k+1;
end
m=length(nodes);
for i=1:m
    for j=1:k-1
        if nodes(i,1)==initial_node_list(j,1)
            nodes(i,3)=initial_node_list(j,2);
        end
        if nodes(i,2)==initial_node_list(j,1)
            nodes(i,4)=initial_node_list(j,2);
        end
    end
end
% nodes
% edges
%% Clustering of nodes based on connectivity
%graphconncomp(nodes)
l1=nodes(:,3).';
l2=nodes(:,4).';
G = graph(l1, l2);
plot(G)
bins=conncomp(G);
a4 = unique(bins);

out2 = [a4,histc(bins(1,:),a4)];
out3=out2(1,length(a4)+1:length(out2));
comp_number=unique(bins);
comp_node=[];
comp_edge={};
g=1;
f=1;
initial_node_list(:,3)=bins.';
for i=1:length(nodes)
    for j=1:length(initial_node_list)
        if initial_node_list(j,1)==nodes(i,1)
            nodes(i,5)=initial_node_list(j,3);
        end
        if initial_node_list(j,1)==nodes(i,2)
            nodes(i,6)=initial_node_list(j,3);
        end
    end
end
 %nodes
 g=1;
 k=1;
 for i=1:length(comp_number)
     f=1;
     for j=1:length(nodes)
         if comp_number(1,i) == nodes(j,5)
             comp_node(f,g)=nodes(j,1);
             comp_node(f,g+1)=nodes(j,2);
             comp_edge(f,k)=edges(j,1);
             f=f+1;
         end
     end
     g=g+2;
     k=k+1;
 end
 %% Computation
 s=1;
 temp_final=zeros(2,2,2); 
 node_sequence_final=zeros(1,2,2);
 bind_final=zeros(2,2,2);
 
 for z=1:length(comp_number)
     %s
     nodes=[];
     edges={};
     [x1,y1]=size(comp_node);
     length(comp_node)
     for y=1:x1
         %y
%          disp('entered');
         if comp_node(y,s)~=0
             nodes(y,1)=comp_node(y,s);
             nodes(y,2)=comp_node(y,s+1);
             edges(y,1)=comp_edge(y,z);
         end
         if comp_node(y,s)==0
             break;
         end
     end
%      nodes
%      edges
     
s=s+2;
j=1;k=1;
non_bind_edge={};
non_bind=[];
bind_edge={};

[lm,ln]=size(nodes);
for i=1:lm
    %temp=0;
    if strcmp(edges(i,1),'compound')==0
        
        non_bind_edge(j,1)=edges(i,1);
        non_bind(j,1)=nodes(i,1);
        non_bind(j,2)=nodes(i,2);
        j=j+1;
               
    end
end

temp_edge=edges;
temp_node=nodes;
edges={};
nodes=[];
edges=non_bind_edge;
nodes=non_bind;
%nodes
%edges

%% Computation Starts
nodecount=unique(nodes);

c=0;
i=1;
node_sequence=zeros(1,2);
temp=zeros(2,2);

if isempty(nodes)==1
    temp_final(1,1,z)=0;
continue;
 end
%% Initial truthtable creation
node_sequence(1,1)=nodes(1,1);
node_sequence(1,2)=nodes(1,2);

    %length(temp)
    i=1;
    if strcmp(edges(1,1),'activation')==1 || strcmp(edges(1,1),'phosphorylation')==1 || strcmp(edges(1,1),'expression')==1
        temp(1,1)=0; temp(1,2)=0;temp(2,1)=0;temp(2,2)=1;temp(3,1)=1;temp(3,2)=1;
        
    end
    if strcmp(edges(1,1),'inhibition')==1 || strcmp(edges(1,1),'dephosphorylation')==1
        temp(1,1)=0;temp(1,2)=0;temp(2,1)=0; temp(2,2)=1;temp(3,1)=1;temp(3,2)=0;
    end
     if strcmp(edges(1,1),'binding/association')==1 
         %temp(1,1)=1; temp(1,2)=1;
         temp(1,1)=0; temp(1,2)=0;temp(2,1)=0;temp(2,2)=1;temp(3,1)=1;temp(3,2)=0;temp(4,1)=1;temp(4,2)=1;
     end
    if strcmp(edges(1,1),'ubiquitnation') == 1
        temp(1,1)=0;temp(1,2)=0;
        temp(2,1)=1;temp(2,2)=1;
    end
    if strcmp(edges(1,1),'dissociation')==1
        temp(1,1)=0;temp(1,2)=0;
        temp(2,1)=1;temp(2,2)=1;
    end


%% Keep Adding New variables
f=3;
h=2;
%disp('ended');
[km,kn]=size(nodes);
flag_bind=0;
while h<=km
    z
    h
    %temp
    %nodes(h,:)
    [m,n]=size(temp);
        
    flag=0;
    if any(node_sequence(1,:) == nodes(h,1))==1 && any(node_sequence(1,:) == nodes(h,2))==0
            disp('match1')
            nodes(h,:)
            flag=1;
            uncommon=nodes(h,2);
            common=nodes(h,1);
            p=length(node_sequence)+1;
            node_sequence(1,p)=nodes(h,2);
            node_to_be_added=h;
            [m,n]=size(temp);
            temp1=[];
            parfor t=1:m
                temp1(t,:)=temp(t,:);
                
            end
            parfor t=1:m
                
                temp1(m+t,:)=temp(t,:);
            end
            temp=[];
            temp=temp1;
            
            %temp
            f=n+1;
            parfor h1=1:m
               temp(h1,f)=0;
               
            end
            parfor h1=1:m
              
               temp(h1+m,f)=1;
            end
            
     
     elseif any(node_sequence(1,:) == nodes(h,2))==1 && any(node_sequence(1,:) == nodes(h,1))==0
            disp('match2')
            nodes(h,:)
            flag=1;
            uncommon=nodes(h,1);
            common=nodes(h,2);
            node_to_be_added=h;
            p=length(node_sequence)+1;
            node_sequence(1,p)=nodes(h,1);
            [m,n]=size(temp);
            temp1=[];
            parfor t=1:m
                temp1(t,:)=temp(t,:);
                
            end
            parfor t=1:m
                
                temp1(m+t,:)=temp(t,:);
            end
            temp=[];
            temp=temp1;
            
            parfor h1=1:m
               temp(h1,f)=0;
               
            end
            parfor h1=1:m
              
               temp(h1+m,f)=1;
            end
            
        elseif any(node_sequence(1,:) == nodes(h,2))==0 && any(node_sequence(1,:) == nodes(h,1))==0
            disp('no match, new incoming node')
            nodes(h,:)
            [g,h1]=size(temp);
            %g
            %h1
           node_to_be_added=h;
           f1=length(node_sequence);
            node_sequence(1,f1+1)=nodes(h,1);
            c_index=f1+1;
            node_sequence(1,f1+2)=nodes(h,2);
            u_index=f1+2;
             if strcmp(edges(h,1),'activation')==1 || strcmp(edges(h,1),'phosphorylation')==1 || strcmp(edges(h,1),'expression')==1
               temp((g+1):(2*g),1:h1)=temp(1:g,1:h1);
               temp((2*g+1):(3*g),1:h1)=temp(1:g,1:h1);
               temp(1:g,h1+1)=0;temp(1:g,h1+2)=0;
               temp((g+1):(2*g),h1+1)=0;temp((g+1):(2*g),h1+2)=1;
               temp((2*g+1):(3*g),h1+1)=1;temp((2*g+1):(3*g),h1+2)=1;
             end
            if strcmp(edges(h,1),'inhibition')==1 || strcmp(edges(h,1),'dephosphorylation')==1
               temp((g+1):(2*g),1:h1)=temp(1:g,1:h1);
               temp((2*g+1):(3*g),1:h1)=temp(1:g,1:h1);
               temp(1:g,h1+1)=0;temp(1:g,h1+2)=0;
               temp((g+1):(2*g),h1+1)=0;temp((g+1):(2*g),h1+2)=1;
               temp((2*g+1):(3*g),h1+1)=1;temp((2*g+1):(3*g),h1+2)=0;
            end
            if strcmp(edges(h,1),'ubiquitination') == 1
               temp((g+1):(2*g),1:h1)=temp(1:g,1:h1);
               temp(1:g,h1+1)=0;temp(1:g,h1+2)=0;
               temp((g+1):(2*g),h1+1)=1;temp((g+1):(2*g),h1+2)=1;
            end
            if strcmp(edges(h,1),'dissociation') == 1
                temp((g+1):(2*g),1:h1)=temp(1:g,1:h1);
               temp(1:g,h1+1)=0;temp(1:g,h1+2)=0;
               temp((g+1):(2*g),h1+1)=1;temp((g+1):(2*g),h1+2)=1;
            end
            if strcmp(edges(h,1),'binding/association') == 1
               temp((g+1):(2*g),1:h1)=temp(1:g,1:h1);
               temp((2*g+1):(4*g),1:h1)=temp(1:2*g,1:h1);
               %temp((g+1):(2*g),1:h1)=temp(1:g,1:h1);
                  
               temp(1:g,h1+1)=0;temp(1:g,h1+2)=0;
               temp((g+1):2*g,h1+1)=0;temp((g+1):2*g,h1+2)=1;
               temp((2*g+1):3*g,h1+1)=1;temp((2*g+1):3*g,h1+2)=0;
               temp((3*g+1):4*g,h1+1)=1;temp((3*g+1):4*g,h1+2)=1;
               
            end
            flag=0;
            %temp
 
     elseif any(node_sequence(1,:) == nodes(h,2))==1 && any(node_sequence(1,:) == nodes(h,1))==1
        disp('both present');
        nodes(h,:)
        common=nodes(h,1);
        uncommon=nodes(h,2);
        node_to_be_added=h;
        flag=1;
        f=f-1;
        flag_bind=1;
     end

    % get the index of the node_sequence array
    if flag==1
        f=f+1;
        [m,n]=size(nodes);
        for g=1:length(node_sequence)
            if node_sequence(1,g)==common
                c_index=g;
            end
            if node_sequence(1,g)==uncommon
                u_index=g;
            end
        end
     c_index
     u_index
    
    if strcmp(edges(node_to_be_added,1),'activation')==1 || strcmp(edges(node_to_be_added,1),'expression')==1
        v=1;
        disp(edges(node_to_be_added,1));
        [m,n]=size(temp);
        %temp
        c=1;
        temp1=[];
        while v<=m
        if (temp(v,c_index)==0 & temp(v,u_index)==0)||(temp(v,c_index)==0 & temp(v,u_index)==1)||(temp(v,c_index)==1 & temp(v,u_index)==1)
            %c=c+1;
            temp1(c,:)=temp(v,:);
            c=c+1;
        else
            %disp('not going');
        end
            v=v+1;
        end
        temp=temp1;
        temp1=[];
        %temp
        %break
       
    
    elseif strcmp(edges(node_to_be_added,1),'phosphorylation')==1 
        v=1;
        disp('phosphorylation');
        [m,n]=size(temp);
        %temp
        c=1;
        temp1=[];
        while v<=m
        if temp(v,u_index)==temp(v,c_index)
            %c=c+1;
            temp1(c,:)=temp(v,:);
            c=c+1;
            else
            %disp('not going');
        end
            v=v+1;
        end
        temp=temp1;
        temp1=[];
        %temp
        %break
    elseif strcmp(edges(node_to_be_added,1),'inhibition')==1 || strcmp(edges(node_to_be_added,1),'repression')==1 
            v=1;
        disp('inhibition');
        [m,n]=size(temp);
        %temp
        c=1;
        temp1=[];
        while v<=m
        if (temp(v,c_index)==0 & temp(v,u_index)==0)||(temp(v,c_index)==0 & temp(v,u_index)==1)||(temp(v,c_index)==1 & temp(v,u_index)==0)
            %c=c+1;
            temp1(c,:)=temp(v,:);
            c=c+1;
            else
            %disp('not going');
        end
            v=v+1;
        end
        temp=temp1;
        temp1=[];
        
        %temp
        %break
       
    
    elseif strcmp(edges(node_to_be_added,1),'dephosphorylation')==1 
            v=1;
        disp('dephosphorylation');
        [m,n]=size(temp);
        c=1;
        %u_index
        %c_index
        temp1=[];
        while v<=m
        if temp(v,u_index)~=temp(v,c_index)
            %c=c+1;
            temp1(c,:)=temp(v,:);
            c=c+1;
            else
            %disp('not going');
        end
            v=v+1;
        end
        temp=temp1;
        temp1=[];
        %temp
        %break
       
   elseif strcmp(edges(node_to_be_added,1),'ubiquitination') == 1
        v=1;
        %temp
        disp('ubinquitnation');
        [m,n]=size(temp);
        
        temp1=[];
        while v<=m
         x=temp(v,u_index)+temp(v,c_index);
        if temp(v,u_index)==temp(v,c_index)
            
            c=c+1;
            temp1(c,:)=temp(v,:);
            c=c+1;
            else
            %disp('not going');
        end
            v=v+1;
        end
        temp=temp1;
        temp1=[];
        
    elseif strcmp(edges(node_to_be_added,1),'dissociation')==1
        v=1;
        disp('dissociation');
       [m,n]=size(temp);
        c=1;
        temp1=[];
        while v<=m
        if temp(v,u_index)==temp(v,c_index)
            %c=c+1;
            temp1(c,:)=temp(v,:);
            c=c+1;
            else
            %disp('not going');
        end
            v=v+1;
        end
        temp=temp1;
        temp1=[];
        %temp
        %break
       end
    end
    
    if flag==0
        f=f+2;
    end
    h=h+1;
 %temp           
end

[m,n]=size(temp);
size_final(z,1)=m;
size_final(z,2)=n;
temp_final(1:m,1:n,z)=temp(1:m,1:n,1);
m=length(node_sequence);
size_final(z,4)=m;
node_sequence_final(1,1:m,z)=node_sequence(1,1:m,1);
flag_bind=0;
 end
disp('out of loop');
[x,y]=size(size_final);
stable_state_count=1;
parfor i=1:x
    if size_final(i,1)>0
     stable_state_count=stable_state_count*size_final(i,1);
    end
end
stable_state_count

final_score=[];
    sum=0.0;
    parfor z1=1:length(out3)
     %sum=sum+(size_final(z1,1)/out3(1,z1));  
     sum=sum+(size_final(z1,1)/(2^size_final(z1,2)));
    end
    fprintf('Final score %d \n', sum/length(out3));

load gong
sound(y,Fs)
%% Mutation states
mute_temp=[];
mute_temp2=[];
s=1;
s1=1;
cyc(1:50)=0;
[x,y]=size(size_final);
for i=1:x
    row=size_final(i,1);
    column=size_final(i,2);
    for a1=1:row-1
        for a2=a1+1:row
            count=0;
            for a3=1:column
                 if temp_final(a1,a3,i)~=temp_final(a2,a3,i)
                    count=count+1; 
                    pos=a3;
                 end
                 if count>1
                        break
                 end
            end
            if count == 1
                mute_temp(s,1)=a1;
                mute_temp(s,2)=a2;
                mute_temp(s,3)=i;
                mute_temp(s,4)=pos;
                
                mute_temp2(s1,1)=a1;
                mute_temp2(s1,2)=a2;
                mute_temp2(s1,3)=i;
                mute_temp2(s1,4)=pos;
                
                s=s+1;
                s1=s1+1;
            end
        end
    end
    k3=size(mute_temp2);
    if isempty(mute_temp2)==0 && k3(1)>2
            
            %i
            %length(mute_temp2)
            
            f=mute_temp2(:,1:2);
            %p=cycleCountBacktrack('edgeList', f, 20);
            parfor i4=1:length(p)
                cyc(i4)=cyc(i4);%+p(i4);
            end
            p=[];
            mute_temp2=[];
            s1=1;
    else
            fprintf('%dth component cannot support mutation \n',i)
    end
    
end
parfor t=1:length(cyc)
    fprintf('Cycles of size %d before perturbation %d \n', t+2,cyc(t));
end
if isempty(mute_temp)==0
    
    f=unique(mute_temp(1:s-1,3));
    
    hold off
    for i=1:length(f)
        x=[];
        k=1;
        pp=size(mute_temp);
        for j=1:pp(1)
            if mute_temp(j,3)==f(i)
                x(k,1)=mute_temp(j,1);
                x(k,2)=mute_temp(j,2);
                k=k+1;
            end
        end
        
        a = unique(x);
        out = [a,histc(x(:),a)];
       
    end
    x=mute_temp(1:s-1,1:2);
    a = unique(x);
    out = [a,histc(x(:),a)];
    a=[];
    x=[];
    x=mute_temp(1:s-1,4);
    a = unique(mute_temp(1:s-1,4));
    out1 = [a,histc(x(:),a)];
else
    disp('No mutation!')
    
end
if isempty(mute_temp)==0
a=[];
a=unique(mute_temp(:,3));
p=length(a);
p1=ceil(p/3);
for i=1:length(a)
    c=1;
    l1=[];
    l2=[];
    pp=size(mute_temp);
    for j=1:pp(1)
    %for j=1:length(mute_temp)
       
       if mute_temp(j,3)==a(i)
           l1(c)=mute_temp(j,1);
           l2(c)=mute_temp(j,2);
           c=c+1;
       end
    end
   
    G1=[];
    tt=strcat('Subplot for ',num2str(a(i)));
    G1 = graph(l1.', l2.');
    subplot(p1,3,i);
    plot(G1);
    title(tt);
end
suptitle('Before perturbation');
end

    %------------------------------------------------------------------------------------
    disp('---------------Perturbation-------------------');
    %dynamic change in network
    gene_num_list=unique(node_sequence_final);
    gene_name_list={};
    for i=2:length(gene_num_list)
        for j=1:length(label_num)
            if gene_num_list(i)==label_num(j)
                gene_name_list(end+1) = label_name(j);
                
            end
        end
    end
    interactions={'activation','inhibition', 'binding/association', 'ubiquitination', 'expression' };
    %--------------------
    [indx,tf] = listdlg('ListString',gene_name_list, 'SelectionMode','single', 'Name', 'Please select the node you want to perturb', 'ListSize',[550,350]);
    [indx1,tf1] = listdlg('ListString',interactions, 'SelectionMode','single', 'Name', 'Please select the type of interaction between the node to be perturbed and the new node', 'ListSize',[550,350]);
    prompt = {'Enter the new incoming node'};
    title1 = 'Input';
    dims =  [1 50];
    definput = {'555'};
    answer = inputdlg(prompt,title1,dims,definput);
    
    x= gene_num_list(indx+1,1);%node target
    y= str2num(cell2mat(answer(1)));%new node
    x
    y
    %--------------------
    score_list=[];
 
    x= gene_num_list(indx+1,1);
    y=555;
    z=interactions(1,indx1); %interaction type
    cyc1(1:50)=0;
    initial_tempfinal=temp_final;
    [a,b,c]=size(node_sequence_final);
    for i2=1:a
        for j=1:b
            for k=1:c
                if node_sequence_final(i2,j,k)==x
                    i1=i2;j1=j;k1=k;
                end
            end
        end
    end
    i=i1;j=j1;k=k1;
    m=size_final(k,1);
    n=size_final(k,2);
    %j'th column in kth component is the target gene
    temp=temp_final(1:m,1:n,k);
    %temp
    for i=1:m
        temp((m+1):(2*m),1:n)=temp(1:m,1:n);
        temp(1:m,n+1)=0;
        temp((m+1):(2*m),n+1)=1;
    end
    temp1=[];
    c=1;
    flag=0;
    %temp
    if isempty(temp)==0
     for i2=1:m
        if strcmp(char(z),'activation')==1 || strcmp(char(z),'ubiquitination')==1 || strcmp(char(z),'phosphorylation')==1 || strcmp(char(z),'dissociation')==1 || strcmp(char(z),'expression')==1
            flag=1;
            if temp(i2,j) == temp(i2,n+1)
                %keep row
                temp1(c,:)=temp(i2,:);
                c=c+1;
            end
        end
        if strcmp(char(z),'inhibition')==1 || strcmp(char(z),'repression')==1 || strcmp(char(z),'dephosphorylation')==1
            
            flag=1;
            if temp(i2,j) ~= temp(i2,n+1)
                %keep row
                temp1(c,:)=temp(i2,:);
                c=c+1;
            end
        end
    end
    if flag==1
        temp=[];
        temp=temp1;
        temp_final1=temp_final;
        temp_final1(1:m,1:n,k)=0;
        [m,n]=size(temp);
        
        temp_final1(1:m,1:n,k)=temp(1:m,1:n);
        size_final1=size_final;
        size_final1(k,1)=m;
        size_final1(k,2)=n;
    else
        %temp
        temp_final1=temp_final;
        temp_final1(1:2*m,1:n+1,k)=temp(1:2*m,1:n+1);
        size_final1=size_final;
        size_final1(k,1)=2*m;
        size_final1(k,2)=n+1;
    end
    
    node_sequence_final(1,n,k)=y;
    else
        size_final1=size_final;
        temp_final1=temp_final;
    end
    %change in mutation
    disp('mutation started');
    mute_temp1=[];
    mute_temp2=[];
    s=1;
    s1=1;
    [x,y]=size(size_final1);
    for i=1:x
        row=size_final1(i,1);
        column=size_final1(i,2);
        f=[];
        for a1=1:row-1
            for a2=a1+1:row
                count=0;
                for a3=1:column
                    if temp_final1(a1,a3,i)~=temp_final1(a2,a3,i)
                        count=count+1;
                        pos=a3;
                    end
                    if count>1
                        break
                    end
                end
                if count == 1
                    mute_temp1(s,1)=a1;
                    mute_temp1(s,2)=a2;
                    mute_temp1(s,3)=i;
                    mute_temp1(s,4)=pos;
                    
                    mute_temp2(s1,1)=a1;
                    mute_temp2(s1,2)=a2;
                    mute_temp2(s1,3)=i;
                    mute_temp2(s1,4)=pos;
                    
                    s=s+1;
                    s1=s1+1;
                end
            end
            
        end
        k3=size(mute_temp2);
        if isempty(mute_temp2)==0 && k3(1)>2
            
            i
            f=mute_temp2(:,1:2);
            %p=cycleCountBacktrack('edgeList', f, 20);
            for i4=1:length(p)
               cyc1(i4)=cyc1(i4);%+p(i4);
            end
            p=[];
            mute_temp2=[];
            s1=1;
        else
            fprintf('After perturbation, %dth component cannot support mutation \n',i)
            
        end
    end
    
    parfor t=1:length(cyc1)
    fprintf('Cycles of size %d after perturbation %d \n', t+2,cyc1(t));
    end
    final_score=[];
    sum=0;
    parfor z1=1:length(out3)
     sum=sum+(size_final1(z1,1)/(2^size_final1(z1,2)));
    end
    fprintf('Final score %d \n', sum/length(out3));
    score_list=[score_list,sum/length(out3)];
    %sum/length(out3)
 figure
a=[];
if isempty(mute_temp1)==0
a=unique(mute_temp1(:,3));
p=length(a);
p1=ceil(p/3);
for i=1:length(a)
    c=1;
    l1=[];
    l2=[];
    pp=size(mute_temp1);
    for j=1:pp(1)
       
       if mute_temp1(j,3)==a(i)
           l1(c)=mute_temp1(j,1);
           l2(c)=mute_temp1(j,2);
           c=c+1;
       end
    end
    
    G1=[];
    G1 = graph(l1.', l2.');
    tt=strcat('Subplot for ',num2str(a(i)));
    subplot(p1,3,i)
    plot(G1)
    title(tt)
end
end

suptitle('After perturbation');
min(score_list)
