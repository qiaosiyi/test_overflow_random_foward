//	 ##     ##   #   #  
//	#  #   ##     # #   
//	#  #     ##    #    
//	 ## #   ##     #    2016.4.18
//
//  ###################
//###################

//read me : how to use :
//#./main
//输出：Num_of_entry.txt 
//每单个数据包随机转发，统计需要增多多少个流表项。
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

// 邻接表的节点
struct AdjListNode {
	int dest;
	int weight;
	struct AdjListNode* next;
};

// 邻接表 结构体
struct AdjList {
	struct AdjListNode *head;  // 指向头节点
};

// 图结构体，V为顶点个数。array为所有的邻接表
struct Graph {
	int V;
	int *sum;
	struct AdjList* array;
};

//创建邻l接表的节点
struct AdjListNode* newAdjListNode(int dest, int weight) {
	struct AdjListNode* newNode = (struct AdjListNode*) malloc(
			sizeof(struct AdjListNode));
	newNode->dest = dest;
	newNode->weight = weight;
	newNode->next = NULL;
	return newNode;
}

//创建一个图，包含V的顶点
struct Graph* createGraph(int V) {
	struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
	graph->V = V;

	graph->sum = (int *)malloc(V * sizeof(int));

	graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));
	int i = 0;
	for (i; i < V; ++i){
		graph->array[i].head = NULL;
		graph->sum[i]=0;
	}


	return graph;
}
void free_graph(struct Graph* src_route){
	int V,i,j;
	struct AdjListNode* free_Point;
	struct AdjListNode* temp_Point;
	V = src_route->V;
	free(src_route->sum);

	i = 0;
	for(i;i<V;i++){
		free_Point = src_route->array[i].head;
		while(free_Point != NULL){
			temp_Point = free_Point;
			free(temp_Point);
			free_Point = free_Point->next;
		}
	}
	free(src_route->array);
	free(src_route);
}

//添加路由路径到图数据结构内，到目的地dest的一系列hop节点名称
void add_route_point(struct Graph* src_route, int dest, int hop){
	struct AdjListNode* newNode = newAdjListNode(hop, 0);
	(src_route->sum[dest])++;//记录路径长度
	newNode->next = src_route->array[dest].head;
	src_route->array[dest].head = newNode;
}


// 添加一个边(无向图)
void addEdge(struct Graph* graph, int src, int dest, int weight) {

	struct AdjListNode* newNode = newAdjListNode(dest, weight);
	newNode->next = graph->array[src].head;
	graph->array[src].head = newNode;

	newNode = newAdjListNode(src, weight);
	newNode->next = graph->array[dest].head;
	graph->array[dest].head = newNode;
}

void deletePoint(struct Graph* graph, int src){//将这个点对应的所有的
	// printf("delPoint %d\n",src);					//边长设为 10000
	struct AdjListNode* delPoint;
	struct AdjListNode* delPoint2;
	int neighber;
	delPoint = graph->array[src].head;

	while(delPoint != NULL){
		delPoint->weight = 10000;
		neighber = delPoint->dest;
		// printf("delete the %d\n", neighber);
		delPoint2 = graph->array[neighber].head;
		delPoint = delPoint->next;
		while(delPoint2 != NULL){
			if(delPoint2->dest == src){
				delPoint2->weight = 10000;
				// printf("delete the %d\n\n", delPoint2->dest);
				break;
			}
			delPoint2 = delPoint2->next;
		}
	}
} 
void reAddPoint(struct Graph* graph,int src){
	if(src>=0 && src<(graph->V)){
		// printf("reAddPoint %d\n",src);
		struct AdjListNode* reAddPoint;
		struct AdjListNode* reAddPoint2;
		int neighber;
		reAddPoint = graph->array[src].head;
		while(reAddPoint != NULL){
		reAddPoint->weight = 1;
		neighber = reAddPoint->dest;
		// printf("reAdd the %d\n", neighber);
		reAddPoint2 = graph->array[neighber].head;
		reAddPoint = reAddPoint->next;
		while(reAddPoint2 != NULL){
			if(reAddPoint2->dest == src){
				reAddPoint2->weight = 1;
				// printf("reAdd the %d\n\n", reAddPoint2->dest);
				break;
			}
			reAddPoint2 = reAddPoint2->next;
		}
	}
	}else{
		printf("(%d) point does not exist in graph!!\n", src);
	}
}

// 最小堆节点
struct MinHeapNode {
	int v;  //下标
	int dist; //距离
};

// 最小堆
struct MinHeap {
	int size;
	int capacity;
	int *pos;     // pos[i]表示顶点i所在的下标
	struct MinHeapNode **array;
};

// 创建一个最小堆节点
struct MinHeapNode* newMinHeapNode(int v, int dist) {
	struct MinHeapNode* minHeapNode = (struct MinHeapNode*) malloc(
			sizeof(struct MinHeapNode));
	minHeapNode->v = v;
	minHeapNode->dist = dist;
	return minHeapNode;
}

// A utility function to create a Min Heap
struct MinHeap* createMinHeap(int capacity) {
	struct MinHeap* minHeap = (struct MinHeap*) malloc(sizeof(struct MinHeap));
	minHeap->pos = (int *) malloc(capacity * sizeof(int));
	minHeap->size = 0;
	minHeap->capacity = capacity;
	minHeap->array = (struct MinHeapNode**) malloc(
			capacity * sizeof(struct MinHeapNode*));
	return minHeap;
}

// 交换两个最小堆的节点
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b) {
	struct MinHeapNode* t = *a;
	*a = *b;
	*b = t;
}

//在位置 idx 调整堆
void minHeapify(struct MinHeap* minHeap, int idx) {
	int smallest, left, right;
	smallest = idx;
	left = 2 * idx + 1;
	right = 2 * idx + 2;

	if (left < minHeap->size
			&& minHeap->array[left]->dist < minHeap->array[smallest]->dist)
		smallest = left;

	if (right < minHeap->size
			&& minHeap->array[right]->dist < minHeap->array[smallest]->dist)
		smallest = right;

	if (smallest != idx) {
		// 需要交换的节点
		struct MinHeapNode *smallestNode = minHeap->array[smallest];
		struct MinHeapNode *idxNode = minHeap->array[idx];

		//交换下标
		minHeap->pos[smallestNode->v] = idx;
		minHeap->pos[idxNode->v] = smallest;

		//交换节点
		swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

		minHeapify(minHeap, smallest);
	}
}

// 推是否为空
int isEmpty(struct MinHeap* minHeap) {
	return minHeap->size == 0;
}

struct MinHeapNode* extractMin(struct MinHeap* minHeap) {
	if (isEmpty(minHeap))
		return NULL;

	struct MinHeapNode* root = minHeap->array[0];

	struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
	minHeap->array[0] = lastNode;

	// 更新下标
	minHeap->pos[root->v] = minHeap->size - 1;
	minHeap->pos[lastNode->v] = 0;

	// 记得减少堆的大小
	--minHeap->size;
	minHeapify(minHeap, 0);

	return root;
}

// 当节点v的距离更新后(变小了)调整堆
void decreaseKey(struct MinHeap* minHeap, int v, int dist) {
	//获取节点 v 在 堆中的下标
	int i = minHeap->pos[v];

	minHeap->array[i]->dist = dist;

	// 因为是变小了，自下向上调整堆即可。 O(Logn)
	while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {
		minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
		minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
		swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

		i = (i - 1) / 2;
	}
}

// 判断节点v是否在堆中
bool isInMinHeap(struct MinHeap *minHeap, int v) {
	if (minHeap->pos[v] < minHeap->size)
		return true;
	return false;
}

// 打印结果
void printArr(int dist[], int n) {
	printf("Vertex   Distance from Source\n");
	int i = 0;
	for (i; i < n; ++i)
		printf("%d \t\t %d\n", i, dist[i]);
	i = 0;

}

void dijkstra(FILE *ofp,int *distance,struct Graph* graph,
				struct Graph* src_route, int src) {

	int V = graph->V;
	//int dist[V];
	int *dist;
	int parent[V];
	int num_of_entry[V];
	parent[src]=src;//the pre point to src is the src it self
	dist = distance;

	struct MinHeap* minHeap = createMinHeap(V);
	struct MinHeapNode *array_temp[V];
	// 初始化堆包含所有的顶点
	int v = 0;
	for (v; v < V; ++v) {
		if(v == src){
			minHeap->array[src] = newMinHeapNode(src, dist[src]);
			array_temp[v] = minHeap->array[src];
			continue;
		}
		dist[v] = INT_MAX;
		minHeap->array[v] = newMinHeapNode(v, dist[v]);
		array_temp[v] = minHeap->array[v];
		minHeap->pos[v] = v;
	}

	// 把 源点 src 的距离设置为0，第一个取出的点即为源点
	
	minHeap->pos[src] = src;
	dist[src] = 0;
	decreaseKey(minHeap, src, dist[src]);

	minHeap->size = V;
	int flag;
	flag = 0;
	// 这个循环中，minHeap包含的是所有未在SPT中的顶点
	while (!isEmpty(minHeap)) {
		// 取得堆顶节点，即最小距离的顶点
		struct MinHeapNode* minHeapNode = extractMin(minHeap);
		int u = minHeapNode->v;
		if (flag == 0){
			//printf("u=%d \n",u);//打印源点
			flag = 1;
		}
		// 只需要遍历和u相邻的顶点进行更新
		struct AdjListNode* pCrawl = graph->array[u].head;
		
		while (pCrawl != NULL) {
			int v = pCrawl->dest;
			// 松弛操作，更新距离
			if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX
					&& pCrawl->weight + dist[u] < dist[v]) {
				dist[v] = dist[u] + pCrawl->weight;
				parent[v]=u;
				//printf("parent[%d]=%d\n", v,u);
				//距离更新了之后，要调整最小堆
				decreaseKey(minHeap, v, dist[v]);
			}
			pCrawl = pCrawl->next;
		}
	}

	// 打印
	// printArr(dist, V);  
	int i,pre;
	// i = 0;

	// // printf("\n\ndest\tpre\n");
	for(i;i<V;++i){
		// printf("%d\t%d\n",i,parent[i]);
		num_of_entry[i]=0;
	}

	i = 0;
	// fprintf(ofp,"route of dijkstra \n");//打印出来各个点到src的中间
	for(i;i<V;++i){						//路径,目的地是i，同一行内是hop。
		if(i != src){
			// fprintf(ofp,"%d ",i );	
			add_route_point(src_route, i, i);//struct Graph* src_route
			num_of_entry[i]++;				//, int dest, int hop
		}
		pre = parent[i];
		while(pre != src){

			// fprintf(ofp,"%d ", pre);
			add_route_point(src_route, i, pre);
			num_of_entry[pre]++;
			pre =parent[pre];
		}
		
		// fprintf(ofp,"%d\n",src);
		add_route_point(src_route, i, src);
		num_of_entry[src]++;
	}							//打印出来各个点到src的中间路径


	// i=0;
	// fprintf(ofp,"\nnumofentry\n");//打印各个节点的流表数量
	// for(i;i<V;i++){
	// 	fprintf(ofp,"%d\n",num_of_entry[i] );
	// }

	v = 0;
	for (v; v < V; v++) {
		free(array_temp[v]);
	}
	free(minHeap->pos);
	free(minHeap->array);
	free(minHeap);
}

//将from图里的hop无重复地加入到to图内
void add_graph_from_to(int src,struct Graph* from, struct Graph* to){
	struct AdjListNode* from_hop_List;
	struct AdjListNode* to_hop_List;
	int V,i,exist,add_hop;
	V=from->V;
	i = 0;
	exist = 0;//首先图to内的hop假设没有相同的，有相同的（exist=1）就不加入。
	for(i;i<V;i++){//需要复制V-1个路径
		if(i != src){//源点不复制
			from_hop_List = from->array[i].head;
			to_hop_List = to->array[i].head;
			if(to_hop_List = NULL){//如果allhop是空的，简单完全复制，
				while(from_hop_List != NULL){//将from的hops加完为止
					add_route_point(to, i, from_hop_List->dest);
					from_hop_List = from_hop_List->next;
				}
				continue;//加完一层后，可以跳过之后的添加。
			}else{//如果allhop不是空的，判断有没有重复，然后复制，
				while(from_hop_List != NULL){//循环要加入的hop次
					add_hop = from_hop_List->dest;
					to_hop_List = to->array[i].head;
					while(to_hop_List != NULL){//循环判断是否有已经加入的hop
						if(to_hop_List->dest == add_hop){
							exist = 1;
							break;
						}
						to_hop_List = to_hop_List->next;
					}
					if(exist == 0){//如果没有重复的，加入图to的hop里。
						add_route_point(to, i, add_hop);
					}
					exist = 0;
					from_hop_List = from_hop_List->next;
				}
			}
		}	
	}	
}
void print_graph(struct Graph* graph){
	int V,i;
	struct AdjListNode* print_List;
	V = graph->V;
	i = 0;
	for(i;i<V;i++){
		print_List = graph->array[i].head;
		while(print_List != NULL){
			printf("%d ", print_List->dest);
			print_List = print_List->next;
		}
		printf("\t\t\tsum=%d\n", graph->sum[i]);
		printf("\n");
	}
}
void sum0graph(struct Graph* graph){
	int V,i;
	V = graph->V;
	i = 0;
	for(i;i<V;i++){
		graph->sum[i] = 0;
	}
}

// 测试
int main() {
	// 创建图
	// ifp = fopen(argv[1],"r");
	FILE	*ifp, *ofp;
	// ifp = fopen("topofile.txt","r");
	ifp = fopen("topoNo1deg3.txt","r");	
	ofp = fopen("Num_of_entry.txt","w");
	int Max_num,Num_of_src_data,i,j,k,n,Num_of_as,Max_num_link,Max_num_as;
	Max_num_link = 100000;
	Max_num_as = 10000;
	Num_of_src_data = 0;
	Num_of_as = 0;
	i=0;
	int src_data[Max_num_link],as[Max_num_as];
	struct AdjListNode* delPoint;struct AdjListNode* Hop_Point;
	struct AdjListNode* neighbor_List;
	int neighber;
	for(i;i<Max_num_link;i++){
		src_data[i]=0;
	}
	i=0;
	for(i;i<Max_num_as;i++){
		as[i]=0;
	}

	// printf("read the topofile!!\n");
	while((fscanf(ifp,"%d",&n) != EOF)) {//just read the src file.
		src_data[Num_of_src_data] = n;//save the src data.
		as[n]++;
		//printf("%d\n", n);
		Num_of_src_data++;

		//printf("%6d%6d\n",Num_of_src_data,n);
	}
	printf("Num_of_src_data=%d\n",Num_of_src_data);
	i = Max_num_as;
	while(i--){
		if(as[i]>1){
			Num_of_as++;
		}
	}
	printf("Num_of_as=%d\n", Num_of_as);
	int V = Num_of_as;
	struct Graph* graph = createGraph(V);
	i = 0;
	for(i;i<Num_of_src_data/2;i++){
		addEdge(graph,src_data[2*i],src_data[2*i+1],1);//add links to graph
	}
	int distance[V],newdistance[V],sum2_entry[1000],sum_entry[1000],sum_distance[1000],percent[1000],derta[16],cantGetDest=0;//distance存放原始距离
	i=0;									//newdistance存放某个点扩散
	for(i;i<16;i++){	//init				//后的新距离。
		derta[i] = 0;						//derta是新老距离差的个数
	}										//derta[0]表示距离不变的频次
	i=0;									//derta[1]表示距离增加1的频次
	for(i;i<V;i++){		//init
		distance[i] = 0;
		newdistance[i] = 0;
	}
	i=0;									//percent[120]表示路径增加到
	for(i;i<1000;i++){		//init			//原来的120%的频次。
		percent[i] = 0;
		sum_distance[i] = 0;
		sum_entry[i] = 0;
		sum2_entry[i] = 0;
	}
	//struct Graph* src_route = createGraph(V);//存放原始路径。
	int numOfNeighber;
	numOfNeighber = 0;
	i = 0;
	for(i;i<V;i++){ //发送扩散的点，首先计算此点到其他个点的距离。
					//然后，删去此点，计算从邻居到各个点的距离。
					//计数增长了多少的频次。

		struct Graph* src_route = createGraph(V);//存放原始路径。
		struct Graph* all_hop = createGraph(V);//统计到目的地所有的hop内容，单源点。
		struct Graph* all2_hop = createGraph(V);//已经建立过一次路径后新增流表。
		if(i%10==0){
			printf("%d of %d. \n",i ,V-1);//打印处理进度
		}
		dijkstra(ofp,distance,graph,src_route, i);//计算原始距离distance
		j=0;
		// printf("distance\n");
		// for(j;j<V;j++){
		// 	printf("%d\n",distance[j] );
		// }
		deletePoint(graph,i);//在图中删去扩散点，然后计算邻居到各个点的新距离
		neighbor_List = graph->array[i].head;//得到删除点的邻居链表
		
		j=0;//打印路径长度以及hops
		// fprintf(ofp,"\n\nthe distances\n");
		for(j;j<V;j++){
			// fprintf(ofp,"%d\n", src_route->sum[j]);//打印了路径长度
			if(j != i){
				sum_distance[src_route->sum[j]]++;//统计路径长度
			}
			// Hop_Point = src_route->array[j].head;//打印存储在结构内的路径hops
			
			// while(Hop_Point != NULL){     //打印存储在结构内的路径hops

			// 	fprintf(ofp, "%d ", Hop_Point->dest);
			// 	Hop_Point = Hop_Point->next;
			// }
			// fprintf(ofp, "\n");//打印存储在结构内的路径hops
		}
		add_graph_from_to(i,src_route,all2_hop);//原始路径先放到all2_hop里
		sum0graph(all2_hop);//将数量清零
		while(neighbor_List != NULL){//遍历邻居
			struct Graph* neighbor_hop = createGraph(V);//从邻居开始，存放所有的hop。
			numOfNeighber++;
			neighber = neighbor_List->dest;
			dijkstra(ofp,newdistance,graph, neighbor_hop, neighber);//计算扩散后的新距离
			add_graph_from_to(i,neighbor_hop,all_hop);//加入到all_hop里
			add_graph_from_to(i,neighbor_hop,all2_hop);//统计新添加的流表

			neighbor_List = neighbor_List->next;
			free_graph(neighbor_hop);
			if(numOfNeighber>1000){//邻居最多只遍历100个
				break;

			}
		}


		// print_graph(all_hop);
 
		j=0;
		for(j;j<V;j++){
			if(j != i){
				sum_entry[all_hop->sum[j]]++;
				sum2_entry[all2_hop->sum[j]]++;
			}
		}


		// printf("\t%d\n",numOfNeighber );
		numOfNeighber = 0;
		reAddPoint(graph,i);//将删去的点加回来
		free_graph(src_route);
		free_graph(all_hop);
		free_graph(all2_hop);


	}	


	fprintf(ofp, "\n\n总entry数量\t\t频次\n");//打印计算结果，增加百分比
	i = 0;
	for(i;i<1000;i++){
		if(sum_distance[i]>0){
			fprintf(ofp, "\t%d\t\t%d\n",i,sum_distance[i] );			
		}

		sum_distance[999]+=sum_distance[i];
	}
	fprintf(ofp, "总频次 = %d\n",sum_distance[999]/2 );

	fprintf(ofp, "\n\n直接扩散后entry总数量\t\t频次\n");//打印计算结果，增加百分比
	i = 0;
	for(i;i<1000;i++){
		if(sum_entry[i]>0){
			fprintf(ofp, "\t%d\t\t%d\n",i,sum_entry[i] );			
		}

		sum_entry[999]+=sum_entry[i];
	}
	fprintf(ofp, "总频次 = %d\n",sum_entry[999]/2 );

	fprintf(ofp, "\n\n再扩散后entry总数量\t\t频次\n");//打印计算结果，增加百分比
	i = 0;
	for(i;i<1000;i++){
		if(sum2_entry[i]>0){
			fprintf(ofp, "\t%d\t\t%d\n",i,sum2_entry[i] );			
		}

		sum2_entry[999]+=sum2_entry[i];
	}
	fprintf(ofp, "总频次 = %d\n",sum2_entry[999]/2 );



	fclose(ifp);
	fclose(ofp);

	return 0;
}
