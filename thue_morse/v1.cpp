#include<iostream>
#include<bitset>
#include<fstream>
using namespace std;

const unsigned int bytes_to_flush=256;

bool IsEvilNumber(unsigned int i);
unsigned int cal_polinomial(unsigned int x,unsigned int *coefficient,unsigned int deg);

int main(){
    unsigned int init_bit=0,num_bytes=0;
    cout<<"Input the inital bit(0/1): ";
    cin>>init_bit;
    cout<<"Input the length of sequence(byte): ";
    cin>>num_bytes;

    char c0=0,c1=0;
    if(init_bit==0){
        c0=(unsigned int)0b0110'1001;
        c1=(unsigned int)0b1001'0110;
    }else{
        c0=(unsigned int)0b1001'0110;
        c1=(unsigned int)0b0110'1001;
    }

    char buffer[bytes_to_flush];
    unsigned int buffer_size=0;

    ofstream outFile("the_thue_morse_seq.bin",ios::out|ios::binary);
    for(unsigned int i=0;i<num_bytes;i++){
        if(IsEvilNumber(i)){
            buffer[buffer_size++]=c1;
        }else{
            buffer[buffer_size++]=c0;
        }
        if(buffer_size==bytes_to_flush){
            outFile.write(buffer,buffer_size);
            buffer_size=0;
        }
    }
    outFile.write(buffer,buffer_size);
    outFile.close();
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
    }
    return sum;
}
