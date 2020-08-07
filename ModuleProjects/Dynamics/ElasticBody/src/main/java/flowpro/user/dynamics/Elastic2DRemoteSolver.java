package flowpro.user.dynamics;

import flowpro.api.*;
import flowpro.user.dynamics.Elastic2DRemoteSolver.Body;
import java.io.*;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.SocketAddress;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;


public class Elastic2DRemoteSolver implements Dynamics {

	Equation eqn;

	private int nBodies;
	private Body[] bodies;
	public double dt;
	public double t;
	public double zLength;
	public double tKick;

	// dynamic
	boolean dynamicComputation = false;
	protected double lRef = 1;
	protected double pRef = 1;
	protected double rhoRef = 1;
	protected double vRef = 1;
	protected double tRef = 1;
	protected double mRef;
	protected double IRef;
	protected double kRef;
	protected double torRef;
	protected double bRef;
	protected double torbRef;

	//Matlab
	RemoteStructureSolver mc;

	@Override
	public void init(int nBodies, String simulationPath, String meshPath, Equation eqn) throws IOException {
		this.eqn = eqn;
		this.nBodies = nBodies;
		bodies = new Body[nBodies];
		for (int i = 0; i < nBodies; i++) {
			bodies[i] = new Body(i);
		}

		// clear file for results
		FlowProProperties props = new FlowProProperties();
		props.load(new FileInputStream(simulationPath + "parameters.txt"));

		// read boundary points
		try {
			double[][] PXY = Mat.loadDoubleMatrix(meshPath + "vertices.txt"); // mesh vertices coordinates
			int[][] boundaryALE = Mat.loadIntMatrix(meshPath + "boundaryTypeALE.txt");
			for (int k = 0; k < nBodies; k++) {
				int[] aux = new int[PXY.length];
				for (int i = 0; i < boundaryALE.length; i++) {
					if (boundaryALE[i][0] == k + 2) {
						aux[boundaryALE[i][1]] = 1;
						aux[boundaryALE[i][2]] = 1;
					}
				}
				int s = 0;
				for (int i = 0; i < aux.length; i++) {
					s += aux[i];
				}
				bodies[k].nBoundary = s;
				bodies[k].radialBasisCoefficients = null;
				bodies[k].boundaryPointsCoords = new double[bodies[k].nBoundary][2];
				s = 0;
				for (int i = 0; i < aux.length; i++) {
					if (aux[i] == 1) {
						bodies[k].boundaryPointsCoords[s] = PXY[i];
						s++;
					}
				}
			}
		} catch (IOException ioe) {
			System.out.println("Error in boundary points!" + ioe);
		}

		zLength = 1;
		if (props.containsKey("zLength")) {
			zLength = props.getDouble("zLength");
		} else {
			System.out.println("Length of 2D bodies in z coordinate is set to " + zLength);
		}

		tKick = Double.MAX_VALUE;
		if (props.containsKey("tKick")) {
			tKick = props.getDouble("tKick");
		}

		try {
			double[] refValues = eqn.getReferenceValues();
			lRef = refValues[0];
			pRef = refValues[1];
			rhoRef = refValues[2];
			tRef = refValues[4];
		} catch (Exception e) {
			System.out.println("Cannot assign referential values to body dynamics!");
		}
		mRef = rhoRef * lRef * lRef;
		kRef = pRef;
		bRef = pRef * tRef;
		IRef = mRef * lRef * lRef;
		torRef = pRef * lRef * lRef;
		torbRef = pRef * lRef * lRef * tRef;

		dynamicComputation = true;

		// Launching Matlab
		try {
			System.out.println("lauching remote structure solver");
//            mc = new MatlabClient();
			mc = new JsonRemoteStructureSolver();
			mc.init();
		} catch (Exception e) {
			System.out.println("Matlab init error " + e);
		}
	}

	@Override
	public void computeBodyMove(double dt, double t, int newtonIter, FluidForces fluFor) throws IOException {
		this.dt = dt;
		this.t = t;

		double[][] boundaryForces = fluFor.getBoundaryForce();
		if (boundaryForces != null) {
			int[][] boundaryIndexes = fluFor.getFaceIndexes();
			int nIndex = boundaryForces.length;
			int dim = boundaryForces[0].length / 2;
			int[] nbIndex = new int[nBodies];

			// number of bodies faces
			for (int i = 0; i < nIndex; i++) {
				nbIndex[boundaryIndexes[i][2] - 2]++;
			}
			for (int k = 0; k < nBodies; k++) {
				double[][] bodyForce = new double[dim][nbIndex[k]];
				double[][] forceCoord = new double[nbIndex[k]][dim];
				int s = 0;
				for (int i = 0; i < nIndex; i++) {
					if (boundaryIndexes[i][2] == k + 2) {
						for (int j = 0; j < dim; j++) {
							bodyForce[j][i] = boundaryForces[s][j];
							forceCoord[i][j] = boundaryForces[s][j + dim];
						}
						s++;
					}
				}				
//				mc.computeDeformation(k, bodies[k], forceCoord, bodyForce, t, dt);
				double[][] rbfForceCoefficient = computeRadialBasisFunctionCoefficients(forceCoord, bodyForce);
				mc.computeDeformation(newtonIter, k, bodies[k], forceCoord, rbfForceCoefficient, t, dt);

				bodies[k].radialBasisCoefficients = computeRadialBasisFunctionCoefficients(bodies[k].boundaryPointsCoords, bodies[k].bodyDeformation);
			}
		} else { // external force free
			for (int k = 0; k < nBodies; k++) {
				mc.computeDeformation(newtonIter, k, bodies[k], null, null, t, dt);
				bodies[k].radialBasisCoefficients = computeRadialBasisFunctionCoefficients(bodies[k].boundaryPointsCoords, bodies[k].bodyDeformation);
			}
		}
	}

	@Override
	public MeshMove[] getMeshMove() {
		MeshMove[] mshMov = new MeshMove[nBodies];
		for (int k = 0; k < nBodies; k++) {
			mshMov[k] = new MeshMove(new double[]{0, 0}, new double[]{0}, bodies[k].getRadialBasisCoefficients(), bodies[k].getBoundaryPointsCoords());
		}
		return mshMov;
	}

	@Override
	public double[][] getCenter() {
		return new double[2][nBodies];
	}

	@Override
	public void nextTimeLevel() throws IOException {
		mc.nextTimeLevel();
	}

	@Override
	public void savePositionsAndForces() {
	}

	double[][] computeRadialBasisFunctionCoefficients(double[][] boundaryPointsCoords, double[][] b) {
		int nPoints = boundaryPointsCoords.length;
		int dim = boundaryPointsCoords[0].length;
		double[][] rbfCoeff = new double[dim][nPoints];
		double[][] A = new double[nPoints][nPoints];
		for (int i = 0; i < nPoints; i++) {
			for (int j = 0; j < nPoints; j++) {
				A[i][j] = radialBasisFunction(boundaryPointsCoords[i], boundaryPointsCoords[j]);
			}
		}
		for (int d = 0; d < dim; d++) {
			rbfCoeff[d] = Mat.cg(A, b[d], 1e-6, 500);
		}
		return rbfCoeff;
	}

	public double radialBasisFunction(double[] a, double[] b) {
		double norm = 0;
		for (int i = 0; i < a.length; i++) {
			norm += (a[i] - b[i]) * (a[i] - b[i]);
		}
		norm = Math.sqrt(norm);
		return 1 - norm;
	}

	public class Body {

		// deformation
		public int nBoundary;
		public double[][] radialBasisCoefficients;
		public double[][] boundaryPointsCoords;
		public double[][] bodyDeformation;

		Body(int i) {
		}

		public double[][] getBoundaryPointsCoords() {
			return boundaryPointsCoords;
		}

		public double[][] getRadialBasisCoefficients() {
			return radialBasisCoefficients;
		}
	}
}

interface RemoteStructureSolver {

	public void init() throws IOException;

	public void computeDeformation(int newtonInter, int bodyNumber, Body body, double[][] forceCoord, double[][] rbfForceCoefficient, double t, double dt) throws IOException;

	public void nextTimeLevel() throws IOException;

	public void close() throws IOException;
}

class MatlabClient implements RemoteStructureSolver {

	Socket socket = null;
	ObjectOutputStream out;
	ObjectInputStream in;

	MatlabClient() throws IOException, ClassNotFoundException {
		try (ServerSocket listener = new ServerSocket(5767)) {
			socket = listener.accept();
			socket.setTcpNoDelay(true);
			socket.setKeepAlive(true);
			out = new ObjectOutputStream(new BufferedOutputStream(socket.getOutputStream()));
			out.writeObject("test");
			out.flush();
			in = new ObjectInputStream(new BufferedInputStream(socket.getInputStream()));
			String testMsg = (String) in.readObject();
			assert testMsg.equals("OK");
			System.out.println("succesfully connected to Matlab ...");
		}
	}

	@Override
	public void init() throws IOException {
		try {
			out.writeObject("init");
			out.flush();
			Object o = in.readObject();
			System.out.println("received: " + o);
		} catch (ClassNotFoundException ex) {
			throw new IOException(ex);
		}
	}

	@Override
	public void computeDeformation(int newtonInter, int bodyNumber, Body body, double[][] forceCoord, double[][] rbfForceCoefficient, double t, double dt) throws IOException {
		try {
			out.writeObject("def");
			out.flush();
			out.writeInt(bodyNumber);
			out.writeObject(forceCoord);
			out.writeObject(rbfForceCoefficient);
			out.writeDouble(t);
			out.writeDouble(dt);
			out.flush();
			body.boundaryPointsCoords = (double[][]) in.readObject();
			body.bodyDeformation = (double[][]) in.readObject();
			Object o = in.readObject();
			System.out.println("received: " + o);
		} catch (ClassNotFoundException ex) {
			throw new IOException(ex);
		}
	}

	@Override
	public void nextTimeLevel() throws IOException {
		try {
			out.writeObject("next");
			out.flush();
			Object o = in.readObject();
			System.out.println("received: " + o);
		} catch (ClassNotFoundException ex) {
			throw new IOException(ex);
		}
	}

	@Override
	public void close() throws IOException {
		socket.close();
	}
}

class JsonRemoteStructureSolver implements RemoteStructureSolver {

	public static String HOSTNAME = "localhost";
	public static int PORT = 5767;
	public static final int TIME_OUT = 30000;

	private final Socket socket;
	private final DataOutputStream out;
	private final DataInputStream in;

	JsonRemoteStructureSolver() throws IOException {
		socket = new Socket();
		SocketAddress addr = new InetSocketAddress(HOSTNAME, PORT);
		socket.connect(addr, TIME_OUT);

		out = new DataOutputStream(socket.getOutputStream());
		in = new DataInputStream(socket.getInputStream());

		System.out.println("sending test message");
		out.writeUTF("test");
		out.flush();
		String testMsg = in.readUTF();
		System.out.println("test message received");
		assert testMsg.equals("OK");
		System.out.println("succesfully connected to remote structural solver");
	}

	@Override
	public void init() throws IOException {
		try {
			JSONObject json = new JSONObject();
			json.put("tag", "init");
			out.writeUTF(json.toString());
			out.flush();
		} catch (JSONException ex) {
			throw new IOException("error occured during comunication with remote structure solver: " + ex.getMessage());
		}
	}

	public static double[][] json2Matrix(JSONArray jsonMatrix) throws JSONException {
		int nRows = jsonMatrix.length();
		double[][] mat = new double[nRows][];
		for (int i = 0; i < nRows; i++) {
			JSONArray jsonArray = jsonMatrix.getJSONArray(i);
			int n = jsonArray.length();
			mat[i] = new double[n];

			for (int j = 0; j < n; j++) {
				mat[i][j] = jsonArray.getDouble(j);
			}
		}

		return mat;
	}

	@Override
	public void computeDeformation(int innerIter, int bodyNumber, Body body, double[][] forceCoord, double[][] rbfForceCoefficient, double t, double dt) throws IOException {
		try {
			JSONObject json = new JSONObject();
			json.put("tag", "forces");
			json.put("t", t);
			json.put("dt", dt);
			json.put("innerIter", innerIter);
			json.put("bodyIndex", bodyNumber);
			json.put("forceCoord", new JSONArray(forceCoord));
			json.put("forceRbfCoef", new JSONArray(rbfForceCoefficient));

			out.writeUTF(json.toString());
			out.flush();

			json = new JSONObject(in.readUTF());

			String tag = "deformation";
			if (!json.get("tag").equals(tag)) {
				throw new IOException("message with tag " + tag
						+ " was expected, instead got message with tag " + json.get("tag"));
			}

			body.boundaryPointsCoords = json2Matrix(json.getJSONArray("coordinates"));
			body.bodyDeformation = json2Matrix(json.getJSONArray("displacement"));

		} catch (JSONException ex) {
			throw new IOException("error occured during comunication with remote structure solver: " + ex.getMessage());
		}
	}

	@Override
	public void nextTimeLevel() throws IOException {
		try {
			JSONObject json = new JSONObject();
			json.put("tag", "next");
			out.writeUTF(json.toString());
			out.flush();
		} catch (JSONException ex) {
			throw new IOException("error occured during comunication with remote structure solver: " + ex.getMessage());
		}
	}

	@Override
	public void close() throws IOException {
		try {
			JSONObject json = new JSONObject();
			json.put("tag", "close");
			out.writeUTF(json.toString());
			out.flush();
		} catch (JSONException ex) {
			throw new IOException("error occured during comunication with remote structure solver: " + ex.getMessage());
		}

		socket.close();
	}
}
