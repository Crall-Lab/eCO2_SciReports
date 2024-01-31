int inPin = 10; //declare input pin
int outPin = 7; // declare output pin to control relay
int val = 0;

// the setup function runs once when you press reset or power the board
void setup() {
  // initialize digital pin LED_BUILTIN as an output.
  pinMode(outPin, OUTPUT);
  pinMode(inPin, INPUT);
  Serial.begin(9600);
}

// the loop function runs over and over again forever
void loop() {

  val = digitalRead(inPin);
  if(val == HIGH){
    Serial.println("CO2 normal");
    digitalWrite(outPin, LOW);
  }
    else{
    Serial.println("CO2 low");
    digitalWrite(outPin, HIGH);
    delay(5000);
    digitalWrite(outPin, LOW);

  }
  delay(25000);
}
