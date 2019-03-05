/*******************************************************************************
* Copyright 2014-2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

#include <assert.h>
#include <stdlib.h>

#include "mpi_hpcg_api.hpp"
#include "mpi.h"

/* Comment the next define out to disable forced cancellation of outstanding
 * requests in HPCG_MPI_Request_free() */
#define _HPCG_MPI_REQUEST_FREE_CANCEL_OUTSTANDING_REQUESTS

#ifdef _OPENMPI

/* OpenMPI (and maybe other MPIs) needs additional marshalling of the values.
 * We create dummy values here and do the marshalling in the HPCG_MPI2MPI_*()
 * functions. */

enum {
	_MPI_COMM_NULL = 2,
	_MPI_COMM_WORLD,
	_MPI_BYTE,
	_MPI_DOUBLE,
	_MPI_DOUBLE_INT,
	_MPI_DOUBLE_PRECISION,
	_MPI_FLOAT,
	_MPI_LONG_LONG_INT,
	_MPI_INT,
	_MPI_MAX,
	_MPI_MAXLOC,
	_MPI_MIN,
	_MPI_SUM,
	_MPI_CHAR,
};

#define HPCG_MPI2MPI_CONST(name) case (- _MPI_##name) : return MPI_##name
#define WRAP_MPI_CONST(name) (- _ ## name)

#else

#define WRAP_MPI_CONST(name) name

#endif

int HPCG_MPI_ANY_SOURCE = MPI_ANY_SOURCE;
int HPCG_MPI_MAX_PROCESSOR_NAME = MPI_MAX_PROCESSOR_NAME;

int HPCG_MPI_COMM_NULL = WRAP_MPI_CONST(MPI_COMM_NULL);
int HPCG_MPI_COMM_WORLD = WRAP_MPI_CONST(MPI_COMM_WORLD);

int HPCG_MPI_CHAR = WRAP_MPI_CONST(MPI_CHAR);
int HPCG_MPI_BYTE = WRAP_MPI_CONST(MPI_BYTE);
int HPCG_MPI_DOUBLE = WRAP_MPI_CONST(MPI_DOUBLE);
int HPCG_MPI_DOUBLE_INT = WRAP_MPI_CONST(MPI_DOUBLE_INT);
int HPCG_MPI_DOUBLE_PRECISION = WRAP_MPI_CONST(MPI_DOUBLE_PRECISION);
int HPCG_MPI_FLOAT = WRAP_MPI_CONST(MPI_FLOAT);
int HPCG_MPI_LONG_LONG_INT = WRAP_MPI_CONST(MPI_LONG_LONG_INT);
int HPCG_MPI_SCALAPACK_INT = WRAP_MPI_CONST(MPI_INT);

int HPCG_MPI_MAX = WRAP_MPI_CONST(MPI_MAX);
int HPCG_MPI_MAXLOC = WRAP_MPI_CONST(MPI_MAXLOC);
int HPCG_MPI_MIN = WRAP_MPI_CONST(MPI_MIN);
int HPCG_MPI_SUM = WRAP_MPI_CONST(MPI_SUM);

void *HPCG_MPI_IN_PLACE = (void *)MPI_IN_PLACE;

HPCG_MPI_Request HPCG_MPI_REQUEST_NULL = (HPCG_MPI_Request)MPI_REQUEST_NULL;
HPCG_MPI_Status *HPCG_MPI_STATUSES_IGNORE = (HPCG_MPI_Status *)MPI_STATUSES_IGNORE;
HPCG_MPI_Status *HPCG_MPI_STATUS_IGNORE = (HPCG_MPI_Status *)MPI_STATUS_IGNORE;

int HPCG_MPI_SIMILAR = MPI_SIMILAR;
int HPCG_MPI_SUCCESS = MPI_SUCCESS;
int HPCG_MPI_UNDEFINED = MPI_UNDEFINED;

/* Utils */
static inline void *HPCG_MPI2MPI_Alloc(size_t len)
{
	if (!len) return NULL;
	void *ptr = malloc(len);
	if (ptr == NULL) abort();
	return ptr;
}

static inline void HPCG_MPI2MPI_Free(void *ptr)
{
	free(ptr);
}

#define HPCG_MPI2MPI_Copy_in(dst, src, count, type) do { \
	for (int i = 0; i < count; i++) \
		dst[i] = *((type *)(src + i)); \
} while (0)

#define HPCG_MPI2MPI_Copy_out(dst, src, count, type) do { \
	for (int i = 0; i < count; i++) \
		*((type *)(dst + i)) = src[i]; \
} while (0)

static inline MPI_Request *HPCG_MPI2MPI_Requests_copy_in(int count,
		HPCG_MPI_Request *hpcg_mpi_requests)
{
	MPI_Request *mpi_requests;
	assert(sizeof(MPI_Request) <= sizeof(HPCG_MPI_Request));

	if (sizeof(MPI_Request) == sizeof(HPCG_MPI_Request))
		mpi_requests = (MPI_Request *)hpcg_mpi_requests;
	else {
		mpi_requests =
			(MPI_Request *)HPCG_MPI2MPI_Alloc(sizeof(MPI_Request) * count);
		HPCG_MPI2MPI_Copy_in(mpi_requests,
				hpcg_mpi_requests, count, MPI_Request);
	}

	return mpi_requests;
}

static inline void HPCG_MPI2MPI_Requests_copy_out(int count,
		HPCG_MPI_Request *hpcg_mpi_requests, MPI_Request *mpi_requests)
{
	if (sizeof(MPI_Request) != sizeof(HPCG_MPI_Request)) {
		HPCG_MPI2MPI_Copy_out(hpcg_mpi_requests,
				mpi_requests, count, MPI_Request);
		HPCG_MPI2MPI_Free(mpi_requests);
	}
}

static inline void HPCG_MPI2MPI_Request_copy_in(MPI_Request *mpi_request,
		HPCG_MPI_Request *hpcg_mpi_request)
{
	HPCG_MPI2MPI_Copy_in(mpi_request, hpcg_mpi_request, 1, MPI_Request);
}

static inline void HPCG_MPI2MPI_Request_copy_out(
		HPCG_MPI_Request *hpcg_mpi_request, MPI_Request *mpi_request)
{
	HPCG_MPI2MPI_Copy_out(hpcg_mpi_request, mpi_request, 1, MPI_Request);
}

static inline MPI_Status *HPCG_MPI2MPI_Statuses_copy_in(int count,
		HPCG_MPI_Status *hpcg_mpi_statuses)
{
	MPI_Status *mpi_statuses;
	assert(sizeof(MPI_Status) <= sizeof(HPCG_MPI_Status));

	if (sizeof(MPI_Status) == sizeof(HPCG_MPI_Status))
		mpi_statuses = (MPI_Status *)hpcg_mpi_statuses;
	else if (hpcg_mpi_statuses == HPCG_MPI_STATUSES_IGNORE)
		mpi_statuses = MPI_STATUSES_IGNORE;
	else {
		mpi_statuses =
			(MPI_Status *)HPCG_MPI2MPI_Alloc(sizeof(MPI_Status) * count);
		HPCG_MPI2MPI_Copy_in(mpi_statuses,
				hpcg_mpi_statuses, count, MPI_Status);
	}

	return mpi_statuses;
}

static inline void HPCG_MPI2MPI_Statuses_copy_out(int count,
		HPCG_MPI_Status *hpcg_mpi_statuses, MPI_Status *mpi_statuses)
{
	if (sizeof(MPI_Status) != sizeof(HPCG_MPI_Status)
			&& mpi_statuses != MPI_STATUSES_IGNORE)
	{
		HPCG_MPI2MPI_Copy_out(hpcg_mpi_statuses,
				mpi_statuses, count, MPI_Status);
		HPCG_MPI2MPI_Free(mpi_statuses);
	} else if (mpi_statuses == MPI_STATUSES_IGNORE)
		hpcg_mpi_statuses = HPCG_MPI_STATUSES_IGNORE;
}

static inline void HPCG_MPI2MPI_Status_copy_in(MPI_Status *mpi_status,
		HPCG_MPI_Status *hpcg_mpi_status)
{
	HPCG_MPI2MPI_Copy_in(mpi_status, hpcg_mpi_status, 1, MPI_Status);
}

static inline void HPCG_MPI2MPI_Status_copy_out(
		HPCG_MPI_Status *hpcg_mpi_status, MPI_Status *mpi_status)
{
	HPCG_MPI2MPI_Copy_out(hpcg_mpi_status, mpi_status, 1, MPI_Status);
}

static inline MPI_Comm HPCG_MPI2MPI_Comm(HPCG_MPI_Comm comm)
{
#ifdef _OPENMPI
	switch (comm) {
		HPCG_MPI2MPI_CONST(COMM_WORLD);
		HPCG_MPI2MPI_CONST(COMM_NULL);
	}
	assert(!"Unknown HPCG_MPI_Comm!");
#else
	/* XXX: we might still want to check that the comm is valid. */
	return comm;
#endif
}

static inline MPI_Datatype HPCG_MPI2MPI_Datatype(HPCG_MPI_Datatype type)
{
#ifdef _OPENMPI
	switch (type) {
		HPCG_MPI2MPI_CONST(CHAR);
		HPCG_MPI2MPI_CONST(BYTE);
		HPCG_MPI2MPI_CONST(DOUBLE);
		HPCG_MPI2MPI_CONST(DOUBLE_INT);
		HPCG_MPI2MPI_CONST(DOUBLE_PRECISION);
		HPCG_MPI2MPI_CONST(FLOAT);
		HPCG_MPI2MPI_CONST(INT);
		HPCG_MPI2MPI_CONST(LONG_LONG_INT);
	}
	assert(!"Unknown HPCG_MPI_Datatype!");
	return MPI_DATATYPE_NULL;
#else
	/* XXX: we might still want to check that the type is valid. */
	return type;
#endif
}

static inline MPI_Op HPCG_MPI2MPI_Op(HPCG_MPI_Op op)
{
#ifdef _OPENMPI
	switch (op) {
		HPCG_MPI2MPI_CONST(MAXLOC);
		HPCG_MPI2MPI_CONST(MAX);
		HPCG_MPI2MPI_CONST(MIN);
		HPCG_MPI2MPI_CONST(SUM);
	}
	/* XXX: we might still want to check that the op is valid. */
	assert(!"Unknown HPCG_MPI_Op!");
#else
	return op;
#endif
}

static inline int HPCG_MPI2MPI_Tag(int tag)
{
	/* XXX: we might still want to check that the op is valid and check for
	 * MPI_TAG_ANY. */
	return tag;
}

static inline int HPCG_MPI2MPI_Root(int root)
{
	/* XXX: we might still want to check that the op is valid and check for
	 * MPI_PROC_NULL and MPI_ROOT. */
	return root;
}

/* XXX: we also may need HPCG_MPI2MPI_Send_buf(void *buf) that would check for
 * MPI_IN_PLACE and MPI_BOTTOM. */

static inline int HPCG_MPI2MPI_Source(int source)
{
#ifdef _OPENMPI
	return (source == HPCG_MPI_ANY_SOURCE) ? MPI_ANY_SOURCE : source;
#else
	return source;
#endif
}

/* MPI wrappers */

int HPCG_MPI_Wait(HPCG_MPI_Request *hpcg_mpi_request,
		HPCG_MPI_Status *hpcg_mpi_status)
{
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
	MPI_Status *mpi_status = HPCG_MPI2MPI_Statuses_copy_in(1, hpcg_mpi_status);
	int ierror = MPI_Wait(&mpi_request, mpi_status);
	HPCG_MPI2MPI_Statuses_copy_out(1, hpcg_mpi_status, mpi_status);
	if (mpi_request == MPI_REQUEST_NULL)
		*hpcg_mpi_request = HPCG_MPI_REQUEST_NULL;
	else
		HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
	return ierror;
}

int HPCG_MPI_Request_free(HPCG_MPI_Request *hpcg_mpi_request)
{
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
#ifdef _HPCG_MPI_REQUEST_FREE_CANCEL_OUTSTANDING_REQUESTS
	/* This hack is required for e.g. SGI MPT which hangs in MPI_Finalize() if
	 * there are outstanding recv requests. These requests originate from e.g.
	 * halo exchange code that eagerly posts recv requests for the next
	 * exchange after the current one is complete. Since the next one may or
	 * may not happen, the last eager recv requests remain outstanding forever.
	 *
	 * We cannot MPI_Wait() for such recv requests here because they have no
	 * matching send requests, so we just cancel and carry on.
	 *
	 * The side effect is that ierror handling is borken. */
	int ierror, flag;
	ierror = MPI_Test(&mpi_request, &flag, MPI_STATUS_IGNORE);
	if (ierror != MPI_SUCCESS)
		return ierror;
	if (!flag) {
		ierror = MPI_Cancel(&mpi_request);
		if (ierror != MPI_SUCCESS)
			return ierror;
	}
	ierror = MPI_Request_free(&mpi_request);
#else
	int ierror = MPI_Request_free(&mpi_request);
#endif
	if (mpi_request == MPI_REQUEST_NULL)
		*hpcg_mpi_request = HPCG_MPI_REQUEST_NULL;
	else
		HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
	return ierror;
}

double HPCG_MPI_Wtime()
{
	return MPI_Wtime();
}

int HPCG_MPI_Finalize()
{
	return MPI_Finalize();
}

int HPCG_MPI_Initmpi(int *argc, char ***argv)
{
#ifdef _USE_THREADED_MPI
	int required = MPI_THREAD_MULTIPLE, provided, rc;
	rc = MPI_Init_thread(argc, argv, required, &provided);
	if (provided != required)
		abort();
	return rc;
#else
	return MPI_Init(argc, argv);
#endif
}

int HPCG_MPI_Get_processor_name(char *name, int *resultlen)
{
	return MPI_Get_processor_name(name, resultlen);
}

int HPCG_MPI_Error_string(int errorcode, char *string, int *resultlen)
{
	return MPI_Error_string(errorcode, string, resultlen);
}

int HPCG_MPI_Gather(void *sendbuf, int sendcount, HPCG_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, HPCG_MPI_Datatype recvtype,
                              int root, HPCG_MPI_Comm comm)
{
    return MPI_Gather(sendbuf, sendcount, HPCG_MPI2MPI_Datatype(sendtype), recvbuf, recvcount, HPCG_MPI2MPI_Datatype(recvtype), root, HPCG_MPI2MPI_Comm(comm));
}

int HPCG_MPI_Irecv(void* buf, int count, HPCG_MPI_Datatype datatype,
		int source, int tag, HPCG_MPI_Comm comm,
		HPCG_MPI_Request *hpcg_mpi_request)
{
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
	int ierror = MPI_Irecv(buf, count, HPCG_MPI2MPI_Datatype(datatype), source,
			HPCG_MPI2MPI_Tag(tag), HPCG_MPI2MPI_Comm(comm), &mpi_request);
	HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
	return ierror;
}

int HPCG_MPI_Barrier(HPCG_MPI_Comm comm)
{
	return MPI_Barrier(HPCG_MPI2MPI_Comm(comm));
}

int HPCG_MPI_Comm_rank(HPCG_MPI_Comm comm, int *rank)
{
	return MPI_Comm_rank(HPCG_MPI2MPI_Comm(comm), rank);
}

int HPCG_MPI_Comm_size(HPCG_MPI_Comm comm, int *size)
{
	return MPI_Comm_size(HPCG_MPI2MPI_Comm(comm), size);
}

int HPCG_MPI_Bcast(void *buffer, int count, HPCG_MPI_Datatype datatype,
		int root, HPCG_MPI_Comm comm)
{
	return MPI_Bcast(buffer, count, HPCG_MPI2MPI_Datatype(datatype),
			HPCG_MPI2MPI_Root(root), HPCG_MPI2MPI_Comm(comm));
}

int HPCG_MPI_Send(void *buf, int count, HPCG_MPI_Datatype datatype,
		int dest, int tag, HPCG_MPI_Comm comm)
{
	return MPI_Send(buf, count, HPCG_MPI2MPI_Datatype(datatype),
			dest, HPCG_MPI2MPI_Tag(tag), HPCG_MPI2MPI_Comm(comm));
}

int HPCG_MPI_Startall(int count, HPCG_MPI_Request *hpcg_mpi_requests)
{
	MPI_Request *mpi_requests = HPCG_MPI2MPI_Requests_copy_in(count,
			hpcg_mpi_requests);
	int ierror = MPI_Startall(count, mpi_requests);
	HPCG_MPI2MPI_Requests_copy_out(count, hpcg_mpi_requests, mpi_requests);
	return ierror;
}

int HPCG_MPI_Waitall(int count, HPCG_MPI_Request *hpcg_mpi_requests,
		HPCG_MPI_Status* hpcg_mpi_statuses)
{
	MPI_Request *mpi_requests = HPCG_MPI2MPI_Requests_copy_in(count,
			hpcg_mpi_requests);
	MPI_Status *mpi_statuses = HPCG_MPI2MPI_Statuses_copy_in(count,
			hpcg_mpi_statuses);
	int ierror = MPI_Waitall(count, mpi_requests, mpi_statuses);
	HPCG_MPI2MPI_Statuses_copy_out(count, hpcg_mpi_statuses, mpi_statuses);
	HPCG_MPI2MPI_Requests_copy_out(count, hpcg_mpi_requests, mpi_requests);
	return ierror;
}

int HPCG_MPI_Allreduce(void *sendbuf, void *recvbuf,
		int count, HPCG_MPI_Datatype datatype, HPCG_MPI_Op op,
		HPCG_MPI_Comm comm)
{
	return MPI_Allreduce(sendbuf, recvbuf, count,
			HPCG_MPI2MPI_Datatype(datatype), HPCG_MPI2MPI_Op(op),
			HPCG_MPI2MPI_Comm(comm));
}

int HPCG_MPI_Allgather(void *sendbuf, int sendcount, HPCG_MPI_Datatype sendtype,
		void *recvbuf, int recvcount, HPCG_MPI_Datatype recvtype,
		HPCG_MPI_Comm comm)
{
	/* No need to translate MPI_IN_PLACE because the wrapper value also has a
	 * pointer type */
	return MPI_Allgather(sendbuf, sendcount, HPCG_MPI2MPI_Datatype(sendtype),
			recvbuf, recvcount, HPCG_MPI2MPI_Datatype(recvtype),
			HPCG_MPI2MPI_Comm(comm));
}

int HPCG_MPI_Recv_init(void* buf, int count, HPCG_MPI_Datatype datatype,
		int source, int tag, HPCG_MPI_Comm comm,
		HPCG_MPI_Request *hpcg_mpi_request)
{
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
	int ierror = MPI_Recv_init(buf, count, HPCG_MPI2MPI_Datatype(datatype),
			HPCG_MPI2MPI_Source(source), HPCG_MPI2MPI_Tag(tag),
			HPCG_MPI2MPI_Comm(comm), &mpi_request);
	HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
	return ierror;
}

int HPCG_MPI_Send_init(void* buf, int count, HPCG_MPI_Datatype datatype,
		int dest ,int tag, HPCG_MPI_Comm comm,
		HPCG_MPI_Request *hpcg_mpi_request)
{
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
	int ierror = MPI_Send_init(buf, count, HPCG_MPI2MPI_Datatype(datatype),
			dest, tag, HPCG_MPI2MPI_Comm(comm), &mpi_request);
	HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
	return ierror;
}


int HPCG_MPI_Test(HPCG_MPI_Request *hpcg_mpi_request, int *flag, HPCG_MPI_Status *hpcg_mpi_status)
{
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
	MPI_Status *mpi_status = HPCG_MPI2MPI_Statuses_copy_in(1, hpcg_mpi_status);
	int ierror = MPI_Test(&mpi_request, flag, mpi_status);
	HPCG_MPI2MPI_Statuses_copy_out(1, hpcg_mpi_status, mpi_status);
	HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
	return ierror;
}

int HPCG_MPI_Cancel(HPCG_MPI_Request *hpcg_mpi_request)
{
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
    int ierror = MPI_Cancel(&mpi_request);
	HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
    return ierror;
}

int HPCG_MPI_Abort(HPCG_MPI_Comm comm, int errorcode) {
    return MPI_Abort(HPCG_MPI2MPI_Comm(comm), errorcode);
}


int HPCG_MPI_Testany(int count, HPCG_MPI_Request *hpcg_mpi_requests, int *indx,
               int *flag, HPCG_MPI_Status *hpcg_mpi_statuses) {
    MPI_Request *mpi_requests = HPCG_MPI2MPI_Requests_copy_in(count,
            hpcg_mpi_requests);
    MPI_Status *mpi_statuses = HPCG_MPI2MPI_Statuses_copy_in(count,
            hpcg_mpi_statuses);
    int ierror = MPI_Testany(count, mpi_requests, indx, flag, mpi_statuses);
    HPCG_MPI2MPI_Statuses_copy_out(count, hpcg_mpi_statuses, mpi_statuses);
    HPCG_MPI2MPI_Requests_copy_out(count, hpcg_mpi_requests, mpi_requests);
    return ierror;
}

int HPCG_MPI_Iallreduce(void *sendbuf, void *recvbuf, int count,
        HPCG_MPI_Datatype datatype, HPCG_MPI_Op op, HPCG_MPI_Comm comm,
        HPCG_MPI_Request *hpcg_mpi_request) {
	MPI_Request mpi_request;
	HPCG_MPI2MPI_Request_copy_in(&mpi_request, hpcg_mpi_request);
    int ierror = MPI_Iallreduce(
        sendbuf, recvbuf, count, HPCG_MPI2MPI_Datatype(datatype), HPCG_MPI2MPI_Op(op), HPCG_MPI2MPI_Comm(comm),
        &mpi_request);
	HPCG_MPI2MPI_Request_copy_out(hpcg_mpi_request, &mpi_request);
    return ierror;
}


