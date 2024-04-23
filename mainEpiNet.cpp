#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ranNumbers.h"

/*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- FUNCTIONS =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// Remember pass random generators as reference (&)
void genWSnet(int *edge, int nEdges, int nNodes, int K, float betaWS, Ran &ranUni)
{
	int ii, jj, kk, ee, auxInt, accept;

	// Conect each node with its neighboring nodes
	for (ii=0; ii<nNodes; ii++) for (kk=0; kk<K; kk++)
	{
		jj = ii + kk + 1;
		if (jj >= nNodes) jj -= nNodes;

		auxInt = kk + ii*K;
		edge[auxInt] = jj;
	}

	if (betaWS == 0.0) return;

	// Rewire the edges
	for (ee=0; ee<nEdges; ee++)
	{
		ii = ee/K;
		jj = edge[ee];

		if (betaWS < 1.0)
			if (ranUni.doub() >= betaWS) continue; // Doesn't rewire

		accept = 0;
		while (!accept)
		{
			jj = ranUni.int32()%nNodes;
			if (jj == ii) continue; // avoid self-loop

			// Find if the edge already exists
			for (kk=0; kk<K; kk++)
			{
				auxInt = kk + ii*K;
				if (edge[auxInt] == jj) break;
				auxInt = kk + jj*K;
				if (edge[auxInt] == ii) break;
				if (kk == K-1) accept = 1; // Accept the edge
			}
		}

		edge[ee] = jj; // Rewire with the new node
	}

	return;
}

//---------------------------------------------------------------------------------//
// Remember pass random generators as reference (&)
void genBAnet(int *edge, int nNodes, int K, Ran &ranUni)
{
	int ii, jj, kk, auxInt;

	// Complete graph with the 2K+1 first nodes
	int initNodes = 2*K + 1;
	for (ii=0; ii<initNodes; ii++) for (kk=0; kk<K; kk++)
	{
		jj = ii + kk + 1;
		if (jj >= initNodes) jj -= initNodes;

		auxInt = kk + ii*K;
		edge[auxInt] = jj;
	}

	int *bagIds, lenBagIds;
	lenBagIds = (nNodes - initNodes)*(K+1) + initNodes;
	bagIds = (int*) malloc(lenBagIds*sizeof(int));
	for (ii=0; ii<initNodes; ii++) bagIds[ii] = ii; // Collect the initial nodes in the bag

	// Connect each new node to previous nodes with
	// a probability dependent on their degrees
	int ee, bb, accept, trial;
	int countIds;
	countIds = initNodes;
	bb = initNodes;
	for (ii=initNodes; ii<nNodes; ii++)
	{
		for (kk=0; kk<K; kk++)
		{
			accept = 0;
			while (!accept)
			{
				trial = ranUni.int32()%countIds;

				jj = bagIds[trial];
				if (jj == ii) continue;

				if (kk == 0) accept = 1; // Accept the edge

				// Find if the edge already exists
				for (ee=0; ee<kk; ee++)
				{
					auxInt = ee + ii*K;
					if (edge[auxInt] == jj) break;
					if (ee == kk-1) accept = 1; // Accept the edge
				}
			}

			ee = kk + ii*K;
			edge[ee] = jj; // Save the edge
			bagIds[bb++] = jj; // Collect the index of the node jj
		}

		bagIds[bb++] = ii;
		countIds += K + 1;
	}
	free(bagIds);

	return;
}

//---------------------------------------------------------------------------------//

void epiSimulation(int *newI_vec, short *nodeStatus, short *nodeInfec,
		int *asympTimeNode, int *infecTimeNode,
		int K, float probInfec, float probDevInfec, float probRandomLD,
		int nNodes, int nEdges, int *edge, int initInfec, int maxDays,
		int flagActLD, int ldStart, int ldEnd, int interval, Ran &ranUni)
{
	// Status (SAIR: Susceptible, Asymptomatic, Infected, Removed)
	// 0:S, 1:A, 2:I, 3:R 
	memset(nodeStatus, 0, nNodes*sizeof(short)); // All nodes are susceptibles
	memset(nodeInfec, 0, nNodes*sizeof(short));

	int nInfec, newI, nAsymp;
	int oldI, daysNewI;

	nInfec = 0;
	newI = 0;
	nAsymp = 0;
	oldI = 0;
	daysNewI = 0;

	int nn;
	int auxInt;
	// Choosing random node for Infected status
	for (nn=0; nn<initInfec; nn++)
	{
		auxInt = ranUni.int32()%nNodes;
		while (nodeStatus[auxInt] != 0) auxInt = ranUni.int32()%nNodes;
		nodeStatus[auxInt] = 2; // Infected
		nInfec++;
		newI++;
	}

	int tt, ii, jj, ee;
	int nContagious;
	short iiStatus, jjStatus;
	short flagLockdown, switchLD, count;
	int time, timeLD;

	flagLockdown = 0;
	switchLD = 0;
	count = 0;
	time = 0;
	timeLD = 0;

	for (tt=0; tt<maxDays; tt++)
	{
		newI_vec[tt] += newI; // Store the new infected nodes per day

		nContagious = nAsymp + nInfec;
		if (nContagious == 0) break;

		if (flagActLD) if (!flagLockdown)
                {
                        if (newI > oldI) daysNewI++;
                        else daysNewI = 0;
                        oldI = newI;
                        if (daysNewI > ldStart)
			{
				flagLockdown = 1; // Activate lockdown once
                                switchLD = 1;
			}
		}

                if (flagLockdown == 1)
                {
                	if (interval > 0)
                	{
                        	count++;
                        	if (count > interval)
                        	{
                                       	count = 0;
                                       	if (switchLD) switchLD = 0;
               	                	else switchLD = 1;
                        	}
                        }

                        timeLD++;
                       	if (timeLD > ldEnd)
                        {
                        	flagLockdown = 2;
                        	switchLD = 0;
                        }
                }

		newI = 0;

		// Identifies the suceptible nodes and determine if they will be infected
		for (ee=0; ee<nEdges; ee++)
		{
			// Lockdown resriction
                        if (switchLD) if (ranUni.doub() < probRandomLD) continue;

			ii = ee/K;
			jj = edge[ee];

			iiStatus = nodeStatus[ii];
			jjStatus = nodeStatus[jj];

			if (iiStatus == 0) // Susceptible
			{
				//if (jjStatus == 1) if (ranUni.doub() <= probInfec)
				if (jjStatus == 1 || jjStatus == 2) if (ranUni.doub() <= probInfec)
				{
					nodeInfec[ii] = 1;
				}
			}

			if (jjStatus == 0) // Susceptible
			{
				//if (iiStatus == 1) if (ranUni.doub() <= probInfec)
				if (iiStatus == 1 || iiStatus == 2) if (ranUni.doub() <= probInfec)
				{
					nodeInfec[jj] = 1;
				}
			}
		}

		// Update states
		for (ii=0; ii<nNodes; ii++)
		{
			iiStatus = nodeStatus[ii];

			if (iiStatus == 0) // Suceptible
			{
				if (nodeInfec[ii] == 1)
				{
					nodeStatus[ii] = 1;
					nAsymp++;
				}
				continue;
			}

			if (iiStatus == 1) // Asymptomatic
			{
				asympTimeNode[ii]--;
				if (asympTimeNode[ii] > 0) continue;
				if (ranUni.doub() <= probDevInfec)
				{
					nodeStatus[ii] = 2; // Infected
					nAsymp--;
					nInfec++;
					newI++;
				}
				else
				{
					nodeStatus[ii] = 3; // Removed
					nAsymp--;
				}
				continue;
			}

			if (iiStatus == 2) // Infected
			{
				infecTimeNode[ii]--;
				if (infecTimeNode[ii] > 0) continue;
				nodeStatus[ii] = 3; // Removed
				nInfec--;
				continue;
			}
		}
	}

	return;
}

/*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- MAIN =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

int main()
{
	//===| PARAMETERS |===//
	short err_flag = 0;
	long seed;
	int netModel, aveD, initInfec, ldStart, ldEnd, interval, maxDays,
	    numSims, flagActLD;
	float xNodes, betaWS;
	float probRandomLD, probInfec, probDevInfec;
	char renglon[200];

	// Number of nodes
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%f", &xNodes);

	// Network model (0: Watts-Strogatz 1: Barabasi-Albert)
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &netModel);

	// Average degree of nodes
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &aveD);

	// Probability of rewired (WS)
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%f", &betaWS);

	// Initial infected nodes
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &initInfec);

	// Probability of developing infection
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%f", &probDevInfec);

	// Probability of infecting
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%f", &probInfec);

	// Maximum number of days to simulate
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &maxDays);

	// Number of simulations to average
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &numSims);

	// Seed for random number
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%ld", &seed);

	// Space
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;

	// Lockdown?
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &flagActLD);

	// Probability to restrict en edge during the LD
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%f", &probRandomLD);

	// Start of lockdown
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &ldStart);

	// End of lockdown
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &ldEnd);

	// Lockdown interval
	if (fgets(renglon, sizeof(renglon), stdin) == NULL) err_flag = 1;
	else sscanf(renglon, "%d", &interval);

	if (err_flag)
	{
		printf("Parameter file error (.data)\n");
		exit (1);
	}

	// Initialize random uniform numbers
	Ran ranUni(12345);
	// Asymptomatic time = 3d +- 1d
	Gammadev gammaA(9.0,3.0,seed); // a = (aveTime/stdTime)^2; b = aveTime/stdTime^2
	// Infected time = 10d +- 3d
	Gammadev gammaI(100.0/9.0,10.0/9.0,seed); // a = (aveTime/stdTime)^2; b = aveTime/stdTime^2

	
	char dirFile[100];
	int tt, ss, nn;
	FILE *fNewI;

	int *edge;
	short *nodeStatus, *nodeInfec;	
	int *asympTime, *infecTime;
	int *newI_vec;

	int nNodes = xNodes;
	int K = aveD/2;
	int nEdges = K*nNodes;

	edge = (int*) malloc(nEdges*sizeof(int));
	nodeStatus = (short*) malloc(nNodes*sizeof(short));
	nodeInfec = (short*) malloc(nNodes*sizeof(short));
	asympTime = (int*) malloc(nNodes*sizeof(int));
	infecTime = (int*) malloc(nNodes*sizeof(int));
	newI_vec = (int*) malloc(maxDays*sizeof(int));

	memset(newI_vec, 0, maxDays*sizeof(int)); 

	for (ss=0; ss<numSims; ss++)
	{
		if (netModel == 1) genBAnet(edge, nNodes, K, ranUni); // BA network
		else genWSnet(edge, nEdges, nNodes, K, betaWS, ranUni); // WS network (beataWS = 1.0 --> ER)

		for (nn=0; nn<nNodes; nn++) asympTime[nn] = gammaA.dev();
		for (nn=0; nn<nNodes; nn++) infecTime[nn] = gammaI.dev();

		epiSimulation(newI_vec, nodeStatus, nodeInfec, asympTime, infecTime,
				K, probInfec, probDevInfec, probRandomLD,
                		nNodes, nEdges, edge, initInfec, maxDays,
				flagActLD, ldStart, ldEnd, interval, ranUni);

		//sprintf(dirFile, "dataSim/newI_%d.dat", ss);
		//fNewI = fopen(dirFile, "w");
		//for (tt=0; tt<maxDays; tt++) fprintf(fNewI, "%d\n", newI_vec[tt]);
		//fclose(fNewI);
		//memset(newI_vec, 0, maxDays*sizeof(int)); 
	}

	fNewI = fopen("aveNewI.dat", "w");
	for (tt=0; tt<maxDays; tt++) fprintf(fNewI, "%d\n", newI_vec[tt]/numSims);
	fclose(fNewI);

	free(newI_vec);
	free(edge);
	free(nodeStatus);
	free(nodeInfec);
	free(asympTime);
	free(infecTime);

	return 0;
}
