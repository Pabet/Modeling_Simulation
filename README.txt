LOGGING (Wichtig)
Wenn das Programm aus einem build-Verzeichnis heraus ausgeführt wird, welches
nicht dem Projekt-Verzeichnis entspricht, muss die Log4cxxConfig.xml-Datei in
dieses build-Verzeichnis gelegt werden, da die Logger sonst nicht konfiguriert
werden und z.B. wahllos Speicheradressen in den logs auftauchen

LOGGING auf dem LRZ
- in CMakeLists.txt -> set(LOGGING OFF)
- keine logging librarie wird geladen und logging statements werden ignoriert

Builden aus einem build-Ordner heraus
1) cmake /"Pfad zum Projekt"/MolSym
2) make

Ausführen
- Programm-Aufruf über shell: ./MolSim input.xml
- input.xml ist ein Beispiel, jeder andere name mit .xml funktioniert

Ausführen von Task 3 Aufgabenblatt4:
-der File "input2.xml" enthält Inputdaten für equilibration von Flüssigkeit
-der File "input3.xml" enthält Inputdaten für Sphere;
-während der equilibration phase wird ein txt-file erstellt, dieser file muss als argsv [2] eingegeben werden
-für die zweite Phase mit Sphere muss das Programm neu gestartet werden mit Argumenten input3.xml und txt.file mit Informationen über Partikel
-die Simulation in zweiter Phase startet in Zeit, wenn die ertste Phase sich beendet hat
-video file "task3a4.ogv" illustriert die zweite Phase des Experiments
-video file "Task2small.ogv" illustriert der kleine Experiment von Task 2

Task 5.5
1)                  ./MolSim 5-5-1input.xml //erstellt Checkpoint .txt-Datei
cooling Argon:      ./MolSim coolingArgon.xml .txt
supercooling Argon: ./MolSim supercoolingArgon.xml .txt

XML-Input-Schema
- example.xml (input) und shapes.xsd müssen im build-Ordner (bei IDE meistens Projektordner) liegen

- äußerer Rahmen, Bezug auf shapes.xsd
    <?xml version="1.0"?>
    <shapes xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
            xsi:noNamespaceSchemaLocation="shapes.xsd">
            ...
    </shapes>  (<- am Ende des Dokuments)

- genau ein settings-Block

    <settings id="simulationSettings">
        <delta-t  t="0.25"  />                          //Zeitschritt-Dauer
        <end-time t="0.5"   />                          //Simulationsdauer
        <factor   val="0.0" />
        <domain-size x="20.0" y="30.0" z="20.0" />      //Größe der Domain als 3d-Vektor
        <rcutoff  val="2.0"                     />      //cutoff-Radius
        <x1-boundary-condition val="reflecting"/>       //x-Achse (unten)
        <x2-boundary-condition val="outflow"   />       //x-Achse (oben)
        <y1-boundary-condition val="periodic"  />       //y-Achse (unten)
        <y2-boundary-condition val="periodic"  />       //y-Achse (oben)
        <z1-boundary-condition val="periodic"  />       //z-Achse (unten)
        <z2-boundary-condition val="periodic"  />       //z-Achse (oben)
        <initial-temperature val="30.0"   />            //Anfangstemperatur
        <n-thermostat t="4.0"             />            //Anzahl der Zeitschritte, nach denen das Thermostat angewandt wird
        <target-temperature val="50.0"    />            //Zieltemperatur
        <temperature-difference val="5.0" />            //Schrittgröße, in der die Temperatur geändert wird
        <gravitation val="9.81"           />            //Schwerebeschleunigung
        <r0 val="2.2"     />                            //average bond length (membrane 5.1)
        <k val="300.0"      />                          //stiffness constant
        <fz-up val="0.8"  />                            //Task 5.1 Fz-Kraft
        <parallelisation-method val = "CoolMuc2" />     //none/CoolMuc2/CoolMuc3  (sets num_threads to 1/28/256;Task 5.3)
        <force-calculation-method  val = "LJ"  />       //LJ/ Smoothed LJ/ gravitation (Task 5.5)
        <rl val="1.9" />                                //rl (Task 5.5)
    </settings>


- (0/1) Checkpoint-Block

    <checkpoint id="checkpointSettings">
    	<write-checkpoint bol="true" />                 //write Checkpoint: true/false
    	<read-checkpoint  bol="true" />                 //read  Checkpoint: true/false
    	<write-checkpoint-time t="5.0" />               //Wann der Checkpoint sein soll
    </checkpoint>


- beliebig viele sphere-Blöcke

	<sphere id="sphere1">
	    <type val="2.0"                           />    //Typ zur Unterscheidung in Paraview (muss einzigartig sein)
	    <epsilon val="5.0"                        />    //Epsilon-Wert (teilchenspezifisch)
    	<sigma val="1.0"                          />    //Sigma-Wert   (teilchenspezifisch)
    	<centre   x="10.0" y="0.0" z="0.0"        />    //Mittelpunkt im Koordinatensystem
	    <radius   val="5.0"                       />
	    <velocity x="0.0" y="0.0" z="0.0"         />    //Geschwindigkeit in 3d-Vektor
	    <mass 	  val="1.0" 	                  />
	    <h		  val="1.025"                     />    //Mesh-Width
	</sphere>


- beliebig viele cuboid-Blöcke

	<cuboid id="some_cuboid">
		<type val="3.0"                            />    //Typ zur Unterscheidung in Paraview (muss einzigartig sein)
	    <epsilon val="3.0"                         />    //Epsilon-Wert (teilchenspezifisch)
        <sigma val="2.0"                           />    //Sigma-Wert   (teilchenspezifisch)
		<right-top-point x="10.0" y="0.0" z="0.0"  />    //Oberer-Rechter-Eckpunkt des Cuboids im Koordinatensystem(ideally should be change with +(doman-size)/2)
		<side-lengths    x="7.0"  y="8.0" z="9.0"   />   //Seitenlängen in allen drei Dimensionen
		<velocity        x="6.0"  y="2.0" z="5.0"   />   //Geschwindigkeit als 3d-Vektor
		<mass		   val="4.0"				   />
		<h             val="1.0"                   />    //Mesh-Width
		<membrane bol="false"                      />    //Cuboid ist eine membrane: true/false
	</cuboid>