#pragma once

#ifndef  _MACRODEFINITION_H_
#define  _MACRODEFINITION_H_

#define D3Q19

#ifdef D2Q9
constexpr int _dim{ 2 };
constexpr int _nLatNodes{ 9 };

#elif defined(D3Q19)
constexpr int _dim{ 3 };
constexpr int _nLatNodes{ 19 };

#elif defined(D3Q27)
constexpr int _dim{ 3 };
constexpr int _nLatNodes{ 27 };

#endif 

#endif
