#include <Wire.h>

#define GLOBAL_CONFIG 0x1A

#define ACCEL_CONFIG 0x1C
#define ACCEL_XH 0x3B
#define ACCEL_XL 0x3C
#define ACCEL_YH 0x3D
#define ACCEL_YL 0x3E
#define ACCEL_ZH 0x3F
#define ACCEL_ZL 0x40

#define TEMP_H 0x41
#define TEMP_L 0x42


#define GYRO_CONFIG 0x1B
#define GYRO_XH 0x43
#define GYRO_XL 0x44
#define GYRO_YH 0x45
#define GYRO_YL 0x46
#define GYRO_ZH 0x47
#define GYRO_ZL 0x48




const double g = 9.81;
double prev_time = 0;
double curr_time = 0;
double degreeToRadian = 3.1415926/180.0;
double radianToDegree = 180.0/3.1415926;
// rate gyros work in rad/s. must keep this in mind
double rotMat[9];         // rotation matrix from local to global, the inverse is the transpose
double rotMatInverse[9];  // calculating these numbers is tricky, because intensive calculations might mean a varying time step
double hMatrix[9];        // the disgusting piece of filth we use to calculate the expected euler angle change from the rate gyros

double omega[3];        // from rate_gyros, the real angular velocity, as measured in the body frame
double accelAngles[3];  // angles fetched from the accelerometer, we use to estimate the approx position of the UAV
double deltaAngles[3];  // calculate the change in angles for the euler matrix, calculated purely from the gyros
double eulerAngles[3];  // the thing we use to calculate the rotation matrix, i.e. final angles
double a[3];            // accelaration vector
double omega_bias[3] = { 0.0, 0.0, 0.0 };
double accel_bias[3] = { 0.0, 0.0, 0.0 };
const int MPU = 0x68;

















///////////////////////////  LINEAR ALGEBRA OPERATIONS on VECTORS AND ROTATION MATRICEs /////////////////


void transpose(double* trans, double* mat, size_t nr, size_t nc) {
  // useful for inverting rotation matrices
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      trans[j * nr + i] = mat[i * nc + j];
    }
  }
  return;
}


void rotateVector(double* rotationMatrix, double* vector, size_t r, size_t c) {
  // the code assumes the range of the transformation is the entire input space, or, nullspace = 0
  double temp[r];  // this must be a stack based machine
  double t;
  for (int i = 0; i < r; i++) {
    t = 0.0;
    for (int j = 0; j < c; j++) {
      t += vector[j] * rotationMatrix[i * c + j];
    }
    temp[i] = t;
  }
  for (int i = 0; i < r; i++) {
    vector[i] = temp[i];
  }
  return;
}

void transformVector(double* transformation, double* vector, double* response, size_t r, size_t c) {
  double t;
  for (int i = 0; i < r; i++) {
    t = 0.0;
    for (int j = 0; j < c; j++) t += transformation[i * c + j] * vector[j];
    response[i] = t;
  }
  return;
}

void adjustRotationMatrix() {
  double psi = eulerAngles[2];    // Yaw along z axis or down
  double theta = eulerAngles[1];  // Pitch along y axis or east
  double phi = eulerAngles[0];    // Roll along x axis or north
  rotMat[0] = cos(psi) * cos(theta);
  rotMat[1] = cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi);
  rotMat[2] = sin(psi) * sin(phi) + cos(psi) * cos(phi) * sin(theta);
  rotMat[3] = sin(psi) * cos(theta);
  rotMat[4] = cos(psi) * cos(phi) + sin(psi) * sin(phi) * sin(theta);
  rotMat[5] = sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
  rotMat[6] = sin(theta);
  rotMat[7] = cos(theta) * sin(phi);
  rotMat[8] = cos(theta) * cos(phi);
  transpose(rotMatInverse, rotMat, 3, 3);
}

void adjustHMatrix() {
  double phi = eulerAngles[0];
  double theta = eulerAngles[1];
  double psi = eulerAngles[2];
    // the middle transformation, i.e. theta must NEVER go to 90 degrees, this makes the h_inv singular, and therefore h is undefined.
    // This phenomenon is called gimbal lock.
    // we will side step this by using either the poisson kinematic equations or by using quaternions
  hMatrix[0] = 1;
  hMatrix[1] = tan(theta) * sin(phi);
  hMatrix[2] = tan(theta) * cos(phi);
  hMatrix[3] = 0;
  hMatrix[4] = cos(phi);
  hMatrix[5] = -sin(phi);
  hMatrix[6] = 0;
  hMatrix[7] = sin(phi) / cos(theta);  // singulariy at theta = pi/2 resolved when phi is zero
  hMatrix[8] = cos(phi) / cos(theta);  // singularity at theta = pi/2, resolved only when phi is also pi/2
}

void initializeAllMatrices() {
  adjustRotationMatrix();
  adjustHMatrix();
  // the above values make the matrix identity
}

//////////////////////////////////////////////////////////////// LINEAR ALGEBRA OPERATIONS END //////////////////////////////////////////////////////////////////









//////////////////////////////////////////////////////////////// SENSOR DATA PROCESSING /////////////////////////////////////////////////////////////////////////




void performEulerKinematics() {
  // code to adjust IMU logic
  // for now, we use the H_MATRIX
  // WARNING: THIS CODE ONLY WORKS WHEN THE ROTATION SEQUENCE PSI->THETA->PHI is used!!!
  // let us just use our values from the gyro for now 

  // first we acquire the angles from the accelerometer. This is under the strict assumption the vehicle is not accelarating in any direction other than straight down 

  accelAngles[2] = accelAngles[2]; // this is the yaw rate, the accelerometer cannot measure this
  accelAngles[1] = -atan(a[0]/(sqrt(a[1]*a[1] + a[2]*a[2]))); // this is the pitch, the returned value is in the range of pi to pi radians
  accelAngles[0] = atan(a[1]/(sqrt(a[0]*a[0] + a[2]*a[2]))); // this is the roll

  // the order of rotations in yaw pitch roll for global to body frame
  
  transformVector(hMatrix, omega, deltaAngles, 3, 3);
  curr_time = millis();
  double delta_time = (curr_time - prev_time)/(1000.0);
  prev_time = curr_time;
  for(int i = 0; i < 3; i++) eulerAngles[i] += (deltaAngles[i]*(delta_time));
  adjustHMatrix();
}

void kalmannFilter(){
  // write the code for a 1D Kalmann filter
  // delta theta is obtained in all frames. and is corrected for. I think. 
  
}





void getIMUData() {
  Wire.beginTransmission(MPU);
  Wire.write(ACCEL_XH);
  uint8_t err = Wire.endTransmission(false);
  if (err != 0) {
    Serial.println(" AN ERROR OCCURED ");
    delay(5000);
  }
  Wire.requestFrom(MPU, 14, true);
  if (Wire.available()) {
    int16_t a_x = (Wire.read() << 8) | Wire.read();
    int16_t a_y = (Wire.read() << 8) | Wire.read();
    int16_t a_z = (Wire.read() << 8) | Wire.read();
    a[1] = (a_x / 4096.0) * g;
    a[0] = (a_y / 4096.0) * g;
    a[2] = (a_z / 4096.0) * g;
    int16_t temp = Wire.read()<<8|Wire.read();
    int16_t g_x = (Wire.read() << 8) | Wire.read();  // read the next 6 regs
    int16_t g_y = (Wire.read() << 8) | Wire.read();
    int16_t g_z = (Wire.read() << 8) | Wire.read();
    omega[0] = (g_x / 65.5)*degreeToRadian;  // roll
    omega[1] = (g_y / 65.5)*degreeToRadian;  // pitch
    omega[2] = (g_z / 65.5)*degreeToRadian;  // yaw
  } else Serial.println("HMMMMM");
}


void calibrateIMU() {
  for (int i = 0; i < 1000; i++) {
    getIMUData();
    omega_bias[0] += omega[0];
    omega_bias[1] += omega[1];
    omega_bias[2] += omega[2];
    accel_bias[0] += a[0];
    accel_bias[1] += a[1];
    accel_bias[2] += a[2];
    delay(1);
  }

  omega_bias[0] /= 1000.0;
  omega_bias[1] /= 1000.0;
  omega_bias[2] /= 1000.0;

  accel_bias[0] /= 1000;
  accel_bias[1] /= 1000;
  accel_bias[2] /= 1000;
  accel_bias[2] -= g;

  Serial.println("Calibration complete");
}

void initializeIMU() {

  // configure the IMU for sensitivity and power management
  Wire.beginTransmission(MPU);
  Wire.write(0x6B);
  Wire.write(0);
  Wire.endTransmission(true);

  Wire.beginTransmission(MPU);  //i2c address of mpu6050
  Wire.write(GLOBAL_CONFIG);    // general config register
  Wire.write(0x05);             // write to config register
  // the above line sets the digital lpf(internal to mpu) to cutoff freq 10hz
  Wire.endTransmission();

  ///////////// setup gyros
  Wire.beginTransmission(MPU);
  Wire.write(GYRO_CONFIG);  // access gyro config register
  Wire.write(0x08);         // write to gyro config register
  // the value 0x08 corresponds to 500 deg/s max gyro(full scale range) sens setting
  Wire.endTransmission();

  /////////////setup accel
  Wire.beginTransmission(MPU);
  Wire.write(ACCEL_CONFIG);
  Wire.write(0x10);
  Wire.endTransmission();


  //// calibrate
  calibrateIMU();

  Serial.println("Pernida Parjkyas");
  Serial.println(omega_bias[0]);
  Serial.println(omega_bias[1]);
  Serial.println(omega_bias[2]);

  Serial.println(accel_bias[0]);
  Serial.println(accel_bias[1]);
  Serial.println(accel_bias[2]);

  delay(10000);
}


void acquireIMUData() {
  // write IMU fetch logic
  getIMUData();
  for (int i = 0; i < 3; i++) {
    omega[i] -= omega_bias[i];
    a[i] -= accel_bias[i];
  }
  return;
}


void scanAvailableI2CDev() {
  Serial.println();
  Serial.println("I2C scanner. Scanning ...");
  byte count = 0;

  Wire.begin();
  for (byte i = 1; i < 120; i++) {
    Wire.beginTransmission(i);
    if (Wire.endTransmission() == 0) {
      Serial.print("Found address: ");
      Serial.print(i, DEC);
      Serial.print(" (0x");
      Serial.print(i, HEX);
      Serial.println(")");
      count++;
      delay(1);  // maybe unneeded?
    }            // end of good response
  }              // end of for loop
  Serial.println("Done.");
  Serial.print("Found ");
  Serial.print(count, DEC);
  Serial.println(" device(s).");
  delay(5000);
}

void setup() {
  Serial.begin(57600);
  Wire.begin();
  scanAvailableI2CDev();
  for(int i = 0; i < 3; i++) accelAngles[i] = 0;
  for(int i = 0; i < 3; i++) eulerAngles[i] = 0;
  // initially the orientation of the UAV is locked at level
  initializeAllMatrices();
  initializeIMU();
  delay(5000);
  prev_time = millis();  
  // put your setup code here, to run once:
}

// okay now to fix this shit

// the IMU sends in readings of angular velocity, which we use to determine the next angle

// calculating the eigenvalues of a rotation matrix is easy 
// for one thing we know that there will be one eigen value with the scaling factor of unity 
// meaning we solve for (A - lambdaI)X = 0 to get the null space. This nullspace should be mapped by one vector , which is the axis of rotation
 

void loop() {
  // put your main code here, to run repeatedly:
  // acquireIMUData();
  // rotateVector(rotMat, omega);
  // rotateVector(rotMat, a);
  // adjustRotationMatrix();
  acquireIMUData();
  performEulerKinematics();

  // Serial.println("Accel Values: ");
  // Serial.print(a[0]);
  // Serial.print(' ');
  // Serial.print(a[1]);
  // Serial.print(' ');
  // Serial.println(a[2]);
  Serial.println("Euler Values");
  Serial.print(accelAngles[0]*radianToDegree);
  Serial.print(' ');
  Serial.print(accelAngles[1]*radianToDegree);
  Serial.print(' ');
  Serial.println(accelAngles[2]*radianToDegree);
  delay(1000);
  // Serial.println("Eulerian Values");
  // Serial.print(eulerAngles[0]);
  // Serial.print(' ');
  // Serial.print(eulerAngles[1]);
  // Serial.print(' ');
  // Serial.println(eulerAngles[2]);
}
