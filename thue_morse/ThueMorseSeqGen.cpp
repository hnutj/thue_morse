#include<iostream>
#include<fstream>
#include<bitset>
#include<string>
using namespace std;

const unsigned int bytes_to_flush=256;

bool IsEvilNumber(unsigned int i);
unsigned int cal_polinomial(unsigned int x,unsigned int *coefficient,unsigned int deg);
void gen_byPolinomial(bool init_bit=false,unsigned int length=1000000u,unsigned int *coefficient,unsigned int deg);

int main(){
    unsigned int init_bit=0,num_bytes=0;
    cout<<"Input the inital bit(0/1): ";
    cin>>init_bit;
    cout<<"Input the length of sequence(byte): ";
    cin>>num_bytes;

    unsigned int deg=0,*coefficient=new unsigned int[deg+1];
    cout<<"Input the sampling polinomal degree and coefficients: ";
    cin>>deg;
    for(unsigned int i=0;i<=deg;i++){
        unsigned int c_index=deg-i;
        cout<<"c"<<c_index<<"= ";
        cin>>coefficient[c_index];
    }

    unsigned int num_seqs=1,step=1;
    cout<<"Input the number of sequences: ";
    cin>>num_seqs;
    cout<<"Input the step of coefficient with max deg: ";
    cin>>step;

    for(unsigned int i=0;i<num_seqs;i++){
        gen_byPolinomial(init_bit,num_bytes,coefficient,deg);
        coefficient[deg]+=step;
    }

    return 0;
}

bool IsEvilNumber(unsigned int i)
{
    unsigned int counter=0;
    while(i>0){
        if(i&1u==1u){
            counter++;
        }
        i=i>>1;
    }
    if(counter&1u==1u){
        return true;
    }
    return false;
}

unsigned int cal_polinomial(unsigned int x, unsigned int *coefficient, unsigned int deg)
{
    unsigned int sum=0;
    for(unsigned int i=deg;i>=0;i--){
        sum=x*sum+coefficient[i];
        if(i==0) break;
    }
    return sum;
}

void gen_byPolinomial(bool init_bit, unsigned int length, unsigned int *coefficient, unsigned int deg)
{
    bitset<8> byte("00000000");
    char buffer[bytes_to_flush];

    unsigned int buffer_size=0;
    unsigned int pos=0;

    string fileName="d"+to_string(deg);
    for(unsigned int i=0;i<=deg;i++){
        fileName+="c"+to_string(deg-i)+"_"+to_string(coefficient[deg-i]);
    }
    fileName+=".bin";
    ofstream outFile(fileName,ios::out|ios::binary);
    for(unsigned int i=0;i<length;i++){
        for(int j=0;j<8;j++){
            pos=cal_polinomial(i*8+j,coefficient,deg);
            if(IsEvilNumber(pos)){
                byte[7-j]=true;
            }else{
                byte[7-j]=false;
            }
        }
        buffer[buffer_size++]=(char)byte.to_ulong();
        if(buffer_size==bytes_to_flush){
            outFile.write(buffer,buffer_size);
            buffer_size=0;
        }
        cout<<"Byte____"<<i<<endl;
    }
    outFile.write(buffer,buffer_size);
    outFile.close();
    cout<<"File has been written and successfully closed!"<<endl;
}
