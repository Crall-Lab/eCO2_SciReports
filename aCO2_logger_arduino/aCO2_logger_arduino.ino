// Use this file to log only the temp, relative, humidity, and CO2 from Adafruit SCD30

// Basic demo for readings from Adafruit SCD30
#include <Adafruit_SCD30.h>
#include "RTClib.h"

Adafruit_SCD30  scd30;
const int chipSelect = 10;

//set up real time clock
RTC_PCF8523 rtc;

//Initialize SD card
#include <SPI.h>
#include <SD.h>

//initialize 
//int co2_sensorPin  = A0;
//float co2_val = 0;
//float co2_val_converted = 0;

void setup(void) {
  Serial.begin(115200);
  while (!Serial) delay(10);     // will pause Zero, Leonardo, etc until serial console opens

  Serial.println("Adafruit SCD30 test!");

  // Try to initialize!
  if (!scd30.begin()) {
    Serial.println("Failed to find SCD30 chip");
    while (1) { delay(10); }
  }
  Serial.println("SCD30 Found!");


  // if (!scd30.setMeasurementInterval(10)){
  //   Serial.println("Failed to set measurement interval");
  //   while(1){ delay(10);}
  // }
  Serial.print("Measurement Interval: "); 
  Serial.print(scd30.getMeasurementInterval()); 
  Serial.println(" seconds");


    // see if the card is present and can be initialized:
  if (!SD.begin(chipSelect)) {
    Serial.println("Card failed, or not present");
    // don't do anything more:
    while (1);
  }

if (! rtc.begin()) {
    Serial.println("Couldn't find RTC");
    Serial.flush();
    while (1) delay(10);
  }
  
    rtc.start();

}

void loop() {
  if (scd30.dataReady()){
    Serial.println("Data available!");

    if (!scd30.read()){ Serial.println("Error reading sensor data"); return; }
        File dataFile = SD.open("logger.csv", FILE_WRITE);
          if (dataFile) {
            //co2_val = analogRead(co2_sensorPin);
            //co2_val_converted = (co2_val/1023)*850 + 350;
            //Serial.print(co2_val);
            //Serial.print(",");
            //Serial.print(co2_val_converted);
            //Get timestamp
            DateTime now = rtc.now();
            dataFile.print(now.year(), DEC);
            dataFile.print(',');
            dataFile.print(now.month(), DEC);
            dataFile.print(',');
            dataFile.print(now.day(), DEC);
            dataFile.print(",");
            dataFile.print(now.hour(), DEC);
            dataFile.print(',');
            dataFile.print(now.minute(), DEC);
            dataFile.print(',');
            dataFile.print(now.second(), DEC);
            dataFile.print(',');
            //Serial.print("Temperature: ");
            dataFile.print(scd30.temperature);
            dataFile.print(",");
            //Serial.println(" degrees C");
            
            //Serial.print("Relative Humidity: ");
            dataFile.print(scd30.relative_humidity);
            dataFile.print(",");
            
            //Serial.print("CO2: ");
            dataFile.print(scd30.CO2, 3);
            dataFile.print(",");
            //dataFile.print(co2_val);
            //dataFile.print(",");
            //dataFile.print(co2_val_converted);

            dataFile.println("");
            dataFile.close();
          }
  } else {
    //Serial.println("No data");
  }

  delay(600000);
  //delay(10000);
}
