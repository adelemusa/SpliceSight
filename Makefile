.PHONY: build up down logs restart front back celery

build:
	docker-compose build

up:
	docker-compose up -d

down:
	docker-compose down

logs:
	docker-compose logs -f

restart:
	docker-compose restart

front:
	cd frontend && npm run dev

back:
	cd backend && uvicorn main:app --reload

celery:
	cd backend && celery -A main.celery worker --loglevel=info
