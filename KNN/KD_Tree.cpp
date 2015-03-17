#include <bits/stdc++.h>
using namespace std;
//k维空间找找到离指定点最近的m个点
#define maxn 2005
#define sqr(x) ((x)*(x))
#define DIM (32*32)
#define ls cur<<1
#define rs cur<<1|1

int n,k,q,m,idx;
map <int,int> mp;

struct POINT
{
    int x[DIM];
    int id;
    //重载，用于按维度划分点集
    bool operator < (const POINT &u) const  
    {  
        return x[idx]<u.x[idx];  
    }
}p[maxn],query_p,poi[15];
typedef pair<int,POINT> IP;
priority_queue <IP> resQue;
struct KD_TREE
{
    POINT ptr[maxn<<2];//记录当前节点代表的点
    int ch[maxn<<2];//记录儿子
    void build(int l,int r,int cur=1,int dep=0)
    {
        if(l > r) return;
        ch[cur] = r-l;
        ch[ls] = ch[rs] = -1;
        idx = dep % k;//确定划分的维度
        int mid = (l+r)>>1;

        nth_element(p+l,p+mid,p+r+1);//以mid为界把小于p[mid]的点归到mid左边,其余归到右边
        ptr[cur] = p[mid];
        build(l,mid-1,ls,dep+1);
        build(mid+1,r,rs,dep+1);
    }
    int getdis(POINT a,POINT b)
    {
        int res = 0;
        for(int i=0;i<k;i++)
            res += sqr(a.x[i]-b.x[i]);
        return res;
    }
    void query(POINT p,int m,int cur=1,int dep=0)
    {
        if(ch[cur] == -1) return;
        IP ip(0,ptr[cur]);
        ip.first = getdis(ptr[cur],p);
        int idx = dep%k;
        bool flag = false;
        int x=ls,y=rs;
        if(p.x[idx] >= ptr[cur].x[idx]) swap(x,y);
        if(ch[x] != -1) query(p,m,x,dep+1);
        if(resQue.size() < m)
        {
            resQue.push(ip);
            flag = true;
        }
        else
        {
            if(ip.first < resQue.top().first)
            {
                resQue.pop();
                resQue.push(ip);
            }
            if(sqr(p.x[idx]-ptr[cur].x[idx]) < resQue.top().first) flag = true;
        }
        if(ch[y] != -1 && flag)
            query(p,m,y,dep+1);
    }
}kd;

int stat[20];
char str[2048];

int main()
{
    time_t cur = clock();
    k = DIM;
    m = 3;
    //puts("read sample...");
    FILE* f = freopen("sample.txt","r",stdin);
    int cnt = 0;
    while(~scanf("%s",str))
    {
        int len = strlen(str);
        p[cnt].id = cnt;
        for(int i=0;i<DIM;++i)
            p[cnt].x[i] = str[i]-'0';
        mp[cnt] = str[DIM]-'0';
        cnt++;
    }
    fclose(f);

    //puts("building...");
    kd.build(0,cnt-1);

    //puts("read test and testing...");
    cnt = 0;
    freopen("test.txt","r",stdin);
    int tot = 0;
    int match = 0;
    while(~scanf("%s",str))
    {
        int len = strlen(str);
        //printf("start process test case %d...\n",cnt++);
        tot++;
        for(int i=0;i<DIM;++i)
            query_p.x[i] = str[i]-'0';

        kd.query(query_p,m);
        memset(stat,0,sizeof(stat));
        for(int i=0;i<m;i++)
        {
            stat[mp[resQue.top().second.id]]++;
            resQue.pop();
        }
        int pos = 0;
        int mx = 0;
        for(int i=0;i<=9;++i)
            if(mx < stat[i])
                mx = stat[i] , pos = i;
        if(pos == str[DIM]-'0')
            match ++;
    }

    fclose(f);

    printf("accuracy = %.2f%%\n",(double)match / tot*100);
    printf("time = %dms\n",clock()-cur);
    return 0;
}
