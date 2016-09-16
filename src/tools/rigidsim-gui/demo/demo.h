#ifndef RIGIDSIM_DEMO_H
#   define RIGIDSIM_DEMO_H

class Demo
{
    public:
        virtual int start() = 0;
};

class QKeyEvent;

/*
 * Demo with a OpenGL display
 */
class XDemo
{
    public:
        virtual int  start() = 0;
        virtual int  step()  = 0;

        virtual void draw()  { }
        virtual bool key_pressed(QKeyEvent*) 
        {  return false; }

        virtual ~XDemo(){}
};

#endif
