/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkScaledSOADataArrayTemplate.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkScaledSOADataArrayTemplate_txx
#define vtkScaledSOADataArrayTemplate_txx

#include "vtkScaledSOADataArrayTemplate.h"

#include "vtkArrayIteratorTemplate.h"
#include "vtkBuffer.h"

#include <cassert>

//-----------------------------------------------------------------------------
template <class ValueType>
vtkScaledSOADataArrayTemplate<ValueType>* vtkScaledSOADataArrayTemplate<ValueType>::New()
{
  VTK_STANDARD_NEW_BODY(vtkScaledSOADataArrayTemplate<ValueType>);
}

//-----------------------------------------------------------------------------
template <class ValueType>
vtkScaledSOADataArrayTemplate<ValueType>::vtkScaledSOADataArrayTemplate()
  : AoSCopy(nullptr)
  , Scale(1)
{
}

//-----------------------------------------------------------------------------
template <class ValueType>
vtkScaledSOADataArrayTemplate<ValueType>::~vtkScaledSOADataArrayTemplate()
{
  for (size_t cc = 0; cc < this->Data.size(); ++cc)
  {
    this->Data[cc]->Delete();
  }
  this->Data.clear();
  if (this->AoSCopy)
  {
    this->AoSCopy->Delete();
    this->AoSCopy = nullptr;
  }
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::SetNumberOfComponents(int val)
{
  this->GenericDataArrayType::SetNumberOfComponents(val);
  size_t numComps = static_cast<size_t>(this->GetNumberOfComponents());
  assert(numComps >= 1);
  while (this->Data.size() > numComps)
  {
    this->Data.back()->Delete();
    this->Data.pop_back();
  }
  while (this->Data.size() < numComps)
  {
    this->Data.push_back(vtkBuffer<ValueType>::New());
  }
}

//-----------------------------------------------------------------------------
template <class ValueType>
vtkArrayIterator* vtkScaledSOADataArrayTemplate<ValueType>::NewIterator()
{
  vtkArrayIterator* iter = vtkArrayIteratorTemplate<ValueType>::New();
  iter->Initialize(this);
  return iter;
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::ShallowCopy(vtkDataArray* other)
{
  SelfType* o = SelfType::FastDownCast(other);
  if (o)
  {
    this->Size = o->Size;
    this->MaxId = o->MaxId;
    this->SetName(o->Name);
    this->SetNumberOfComponents(o->NumberOfComponents);
    this->CopyComponentNames(o);
    this->Scale = o->Scale;
    assert(this->Data.size() == o->Data.size());
    for (size_t cc = 0; cc < this->Data.size(); ++cc)
    {
      vtkBuffer<ValueType>* thisBuffer = this->Data[cc];
      vtkBuffer<ValueType>* otherBuffer = o->Data[cc];
      if (thisBuffer != otherBuffer)
      {
        thisBuffer->Delete();
        this->Data[cc] = otherBuffer;
        otherBuffer->Register(nullptr);
      }
    }
    this->DataChanged();
  }
  else
  {
    this->Superclass::ShallowCopy(other);
  }
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::InsertTuples(
  vtkIdType dstStart, vtkIdType n, vtkIdType srcStart, vtkAbstractArray* source)
{
  // First, check for the common case of typeid(source) == typeid(this). This
  // way we don't waste time redoing the other checks in the superclass, and
  // can avoid doing a dispatch for the most common usage of this method.
  SelfType* other = vtkArrayDownCast<SelfType>(source);
  if (!other)
  {
    // Let the superclass handle dispatch/fallback.
    this->Superclass::InsertTuples(dstStart, n, srcStart, source);
    return;
  }

  if (n == 0)
  {
    return;
  }

  int numComps = this->GetNumberOfComponents();
  if (other->GetNumberOfComponents() != numComps)
  {
    vtkErrorMacro("Number of components do not match: Source: "
      << other->GetNumberOfComponents() << " Dest: " << this->GetNumberOfComponents());
    return;
  }

  vtkIdType maxSrcTupleId = srcStart + n - 1;
  vtkIdType maxDstTupleId = dstStart + n - 1;

  if (maxSrcTupleId >= other->GetNumberOfTuples())
  {
    vtkErrorMacro("Source array too small, requested tuple at index "
      << maxSrcTupleId << ", but there are only " << other->GetNumberOfTuples()
      << " tuples in the array.");
    return;
  }

  vtkIdType newSize = (maxDstTupleId + 1) * this->NumberOfComponents;
  if (this->Size < newSize)
  {
    if (!this->Resize(maxDstTupleId + 1))
    {
      vtkErrorMacro("Resize failed.");
      return;
    }
  }

  this->MaxId = std::max(this->MaxId, newSize - 1);

  std::vector<ValueType> vals(numComps);
  for (vtkIdType i = 0; i < n; i++)
  {
    other->GetTypedTuple(i + srcStart, vals.data());
    this->SetTypedTuple(i + dstStart, vals.data()); // will automatically scale data
  }
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::FillTypedComponent(int compIdx, ValueType value)
{
  ValueType* buffer = this->Data[compIdx]->GetBuffer();
  value /= this->Scale;
  std::fill(buffer, buffer + this->GetNumberOfTuples(), value);
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::SetArray(
  int comp, ValueType* array, vtkIdType size, bool updateMaxId, bool save, int deleteMethod)
{
  const int numComps = this->GetNumberOfComponents();
  if (comp >= numComps || comp < 0)
  {
    vtkErrorMacro("Invalid component number '"
      << comp
      << "' specified. "
         "Use `SetNumberOfComponents` first to set the number of components.");
    return;
  }

  this->Data[comp]->SetBuffer(array, size);

  if (deleteMethod == VTK_DATA_ARRAY_DELETE)
  {
    this->Data[comp]->SetFreeFunction(save != 0, ::operator delete[]);
  }
  else if (deleteMethod == VTK_DATA_ARRAY_ALIGNED_FREE)
  {
#ifdef _WIN32
    this->Data[comp]->SetFreeFunction(save != 0, _aligned_free);
#else
    this->Data[comp]->SetFreeFunction(save != 0, free);
#endif
  }
  else if (deleteMethod == VTK_DATA_ARRAY_USER_DEFINED || deleteMethod == VTK_DATA_ARRAY_FREE)
  {
    this->Data[comp]->SetFreeFunction(save != 0, free);
  }

  if (updateMaxId)
  {
    this->Size = numComps * size;
    this->MaxId = this->Size - 1;
  }
  this->DataChanged();
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::SetArrayFreeFunction(void (*callback)(void*))
{
  const int numComps = this->GetNumberOfComponents();
  for (int i = 0; i < numComps; ++i)
  {
    this->SetArrayFreeFunction(i, callback);
  }
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::SetArrayFreeFunction(
  int comp, void (*callback)(void*))
{
  const int numComps = this->GetNumberOfComponents();
  if (comp >= numComps || comp < 0)
  {
    vtkErrorMacro("Invalid component number '"
      << comp
      << "' specified. "
         "Use `SetNumberOfComponents` first to set the number of components.");
    return;
  }
  this->Data[comp]->SetFreeFunction(false, callback);
}

//-----------------------------------------------------------------------------
template <class ValueType>
typename vtkScaledSOADataArrayTemplate<ValueType>::ValueType*
vtkScaledSOADataArrayTemplate<ValueType>::GetComponentArrayPointer(int comp)
{
  const int numComps = this->GetNumberOfComponents();
  if (comp >= numComps || comp < 0)
  {
    vtkErrorMacro("Invalid component number '" << comp << "' specified.");
    return nullptr;
  }

  return this->Data[comp]->GetBuffer();
}

//-----------------------------------------------------------------------------
template <class ValueType>
bool vtkScaledSOADataArrayTemplate<ValueType>::AllocateTuples(vtkIdType numTuples)
{
  for (size_t cc = 0, max = this->Data.size(); cc < max; ++cc)
  {
    if (!this->Data[cc]->Allocate(numTuples))
    {
      return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
template <class ValueType>
bool vtkScaledSOADataArrayTemplate<ValueType>::ReallocateTuples(vtkIdType numTuples)
{
  for (size_t cc = 0, max = this->Data.size(); cc < max; ++cc)
  {
    if (!this->Data[cc]->Reallocate(numTuples))
    {
      return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
template <class ValueType>
void* vtkScaledSOADataArrayTemplate<ValueType>::GetVoidPointer(vtkIdType valueIdx)
{
  // Allow warnings to be silenced:
  const char* silence = getenv("VTK_SILENCE_GET_VOID_POINTER_WARNINGS");
  if (!silence)
  {
    vtkWarningMacro(<< "GetVoidPointer called. This is very expensive for "
                       "non-array-of-structs subclasses, as the scalar array "
                       "must be generated for each call. Using the "
                       "vtkGenericDataArray API with vtkArrayDispatch are "
                       "preferred. Define the environment variable "
                       "VTK_SILENCE_GET_VOID_POINTER_WARNINGS to silence "
                       "this warning. Additionally, for the vtkScaledSOADataArrayTemplate "
                       "class we also set Scale to 1 since we've scaled how "
                       "we're storing the data in memory now. ");
  }

  size_t numValues = this->GetNumberOfValues();

  if (!this->AoSCopy)
  {
    this->AoSCopy = vtkBuffer<ValueType>::New();
  }

  if (!this->AoSCopy->Allocate(static_cast<vtkIdType>(numValues)))
  {
    vtkErrorMacro(<< "Error allocating a buffer of " << numValues << " '"
                  << this->GetDataTypeAsString() << "' elements.");
    return nullptr;
  }

  this->ExportToVoidPointer(static_cast<void*>(this->AoSCopy->GetBuffer()));

  // This is the hacky thing with this class that we now need to set the scale
  // to 1 since we internally are storing the memory in an unscaled manner
  this->Scale = 1.0;

  return static_cast<void*>(this->AoSCopy->GetBuffer() + valueIdx);
}

//-----------------------------------------------------------------------------
template <class ValueType>
void vtkScaledSOADataArrayTemplate<ValueType>::ExportToVoidPointer(void* voidPtr)
{
  vtkIdType numTuples = this->GetNumberOfTuples();
  if (this->NumberOfComponents * numTuples == 0)
  {
    // Nothing to do.
    return;
  }

  if (!voidPtr)
  {
    vtkErrorMacro(<< "Buffer is nullptr.");
    return;
  }

  ValueType* ptr = static_cast<ValueType*>(voidPtr);
  for (vtkIdType t = 0; t < numTuples; ++t)
  {
    for (int c = 0; c < this->NumberOfComponents; ++c)
    {
      *ptr++ = this->Data[c]->GetBuffer()[t] * this->Scale;
    }
  }
}

#endif
