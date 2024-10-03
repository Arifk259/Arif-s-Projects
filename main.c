#include <msp430.h> 
#include <math.h>
#include <stdlib.h>

#define GREEN_LED (BIT0)

/**
 * main.c
 */

/*
 * P1.1 - RXD
 * P1.2 - TXD
 * P1.4 - temperature input
 * P2.3 - IN2
 * P2.4 - IN1
 * P2.6 - ENA
 */

#define NPOINTS 400

void Init(void);
void Init_ADC(void);
void Init_UART(void);
unsigned char getc(void);
void putc(unsigned char c);
void Sample(int n);
void Send(int n);
void setPWM1(float pwm);
void setPWM2(float pwm);


//--------------------------------------------------------
//GlobalVariables
//--------------------------------------------------------

//This is where you store the ADC10MEM readings
volatile unsigned char temps[NPOINTS];


//--------------------------------------------------------
//Main
//--------------------------------------------------------

/*
void main(void)
{
    WDTCTL = WDTPW | WDTHOLD;   // stop watchdog timer
    P1DIR |= BIT6;
    P1SEL |= BIT6;

    while (1){
        setPWM(0.5);
    }
}
*/


void main(void)
{
    // Initialization
    //Stop watchdog timer to prevent time out reset
    WDTCTL = WDTPW | WDTHOLD;

    //Set system clock to calibrated values for 16 MHz
    DCOCTL = CALDCO_16MHZ;
    BCSCTL1 = CALBC1_16MHZ;

    //Set P2.3, P2.4, P2.6 to output
    P2DIR |= BIT3 + BIT4 + BIT6;

    // debugging
    P1DIR |= BIT6;

    // Use table 16 of the datasheet to configure P2SEL and P2SEL2 for P2.6 to have TA0.1 output
    P1SEL |= BIT6;


    //Clear all outputs.
    P1OUT = 0x00;
    P2OUT = 0x00;

    Init_UART();
    Init_ADC();


    //We wait to receive a one-time transmission from MATLAB to begin.
    getc();

    while (1)
    {


        //char isHot = getc(); //Temperature is greater than setting
        char normpwm = getc();
        volatile float pwm = ((float)normpwm - 127)/128;
        //char pwm = getc();
        //abs(pwm));
        //Set TEC to heating mode
        //pwm = 2;
        //if (isHot == 0x30) //pwm < 0)
        if (pwm > 0) //pwm < 0)
        {
            //pwm = pwm/2;
            //Set P2.3 & turn off P2.4
            P2OUT &= ~BIT4;
            //setPWM1(abs(pwm));
            setPWM1(pwm); //Set TEC to heating mode
            //P1OUT |= GREEN_LED;
        }
        else
        { //Set TEC to cooling mode
            //if (isHot == 0x31) //pwm > 0)
            if (pwm < 0) //pwm > 0)
            {
                pwm = pwm*(-1.0); //2.0 + 0.5;
                P2OUT |= BIT4;
                //setPWM2(abs(pwm));
                setPWM2(pwm); //Set TEC to cooling mode
                //P1OUT & ~GREEN_LED;
            }
        }
        Sample(NPOINTS);
        Send(NPOINTS);

    }
}


//--------------------------------------------------------
//Initializing the 10-bit ADCModule
//--------------------------------------------------------
void Init_ADC(void)
{
    //Set up the ADC10CTL1 register to use CONSEQ_2 (Mode 2 - Repeat single channel) and select ADC10 on P1.4.
    ADC10CTL1 |= INCH_4 + CONSEQ_2;

    //Set up the ADC10AE0 register to enable ADC10 on P1.4
    ADC10AE0 |= BIT4;

    //Set up the ADC10CTL0 register and select 4 ADC10CLK cycles sample-and-hold time, select multiple sample and conversion (MSC), and turn on the ADC10.
    ADC10CTL0 |= ADC10SHT_0 + MSC + ADC10ON;

    //Set up the ADC10CTL0 register to Start Conversion and Enable Conversion.
    ADC10CTL0 |= ADC10SC + ENC;
}

//--------------------------------------------------------
//UARTModule
//--------------------------------------------------------

void Init_UART(void)
{
    //Initialize the USCI
    //RXD is on P1.1
    //TXD is on P1.2

    //Configure P1.2 and P1.2 for secondary peripheral function
    P1SEL |= BIT1 + BIT2;
    P1SEL2 |= BIT1 + BIT2;

    /*
    Set the prescalers to hold data with 115200 baud rate.
    Because we have a 16MHz  clock, operating at 115200 baud rate will mean 16MHz/115200 ~ 138.
    The 16-bit value of (UCAxBR0 + UCAxBR1 × 256) forms the prescaler value.
    See table 15-4 from the the user's guide.
    */
    UCA0BR0 = 138; //0x8A;
    UCA0BR1 = 0;

    //Set up the UCA0CTL1 register to use the SMCLK.
    UCA0CTL1 = 0xC0;

    //Set up the UCA0CTL1 to release UART RESET
    UCA0CTL1 &= ~UCSWRST;
}

//--------------------------------------------------------
//Miscellaneous Functions
//--------------------------------------------------------


void Sample(int n)
{
    //Write a for loop and store the ADC10MEM values in the global variable you declared in the beginning.
    //Remember to perform bit-shifting to transfer most of the 10-bit data to your 8-bit variable.
    int i = 0;
    for (i = 0; i < n; i++)
    {
        temps[i] = ADC10MEM>>2;
    }
}


void Send(int n)
{
    //Transmit the values stored in the global variable.
    int i = 0;
    for (i = 0; i < n; i++)
    {
        while(!(IFG2 & UCA0TXIFG));
        //Send character passed in c
        UCA0TXBUF = temps[i];
    }
}

unsigned char getc(void)
{
    //Wait until a character has been received
    while(!(IFG2 & UCA0RXIFG));
    //Fetch the character
    unsigned char c = UCA0RXBUF;
    return c;

}

void setPWM1(float pwm){
    // Determine the value of TACCR0 for a 2 ms pulse
    TACCR0 = 3999;

    // Select an appropriate output mode based on Table 12-2 of the user's guide
    TACCTL1 = OUTMOD_7;

    // Determine the value of TACCR1 for 75% duty cycle
    /*if (pwm > 1)
    {
        pwm = 1;
    }*/
    //int num = floor(7999*pwm);

    TACCR1 = floor(3999 * pwm);

    //Set timerA to TASSEL_2 (SMCLK), and up mode (MC_1: count-to-TACCRO).
    TACTL = TASSEL_2 | MC_1 | ID_2;
}

void setPWM2(float pwm){
    // Determine the value of TACCR0 for a 2 ms pulse
    TACCR0 = 3999;

    // Select an appropriate output mode based on Table 12-2 of the user's guide
    TACCTL1 = OUTMOD_3;

    // Determine the value of TACCR1 for 75% duty cycle
        /*if (pwm > 1)
        {
            pwm = 1;
        }*/
        //int num = floor(7999*pwm);

    TACCR1 = floor(3999 * pwm);

    //Set timerA to TASSEL_2 (SMCLK), and up mode (MC_1: count-to-TACCRO).
    TACTL = TASSEL_2 | MC_1 | ID_2;
}


// Calcuations for PID: This will go in the while loop in main function


